/****************************************************************************************************************************
 * Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Author: Mahamadou Diakite, PhD.
 * Institution: National Institutes of Health (NIH)
 * Lang: C++
 * Date: 10/15/2017
 * Version: 1.0.0
 ****************************************************************************************************************************/

#include "bartgadget.h"
#include <boost/tokenizer.hpp>
#include <sstream>
#include <utility>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <random>
#include <functional>
#include <mutex>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>

#include "bart_api.h"


namespace internal {
     void cleanup(const std::string&);
     class ScopeGuard
     {
     public:
	  ScopeGuard(std::string p) : p_(std::move(p)) {}
	  ~ScopeGuard()
	       {
		    if (is_active_) {
			 cleanup(p_);
		    }
		    deallocate_all_mem_cfl();
	       }

	  void dismiss() { is_active_ = false; }
     private:
	  bool is_active_;
	  const std::string p_;
     };

     
     void cleanup(const std::string &createdFiles)
     {
	  boost::filesystem::remove_all(createdFiles);
     }

     void ltrim(std::string &str)
     {
	  str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](int s) {return !std::isspace(s); }));
     }

     void rtrim(std::string &str)
     {
	  str.erase(std::find_if(str.rbegin(), str.rend(), [](int s) {return !std::isspace(s);}).base(), str.end());
     }

     void trim(std::string &str)
     {
	  ltrim(str);
	  rtrim(str);
     }
	
     std::string get_output_filename(const std::string& bartCommandLine)
     {
	  boost::char_separator<char> sep(" ");
#if 0
	  boost::tokenizer<boost::char_separator<char>,
			   std::string::reverse_iterator> tokens(bartCommandLine.rbegin(),
								 bartCommandLine.rend());
	  auto tok(*tokens.begin());
	  return std::string(tok.rbegin(), tok.rend());
#else 
	  std::vector<std::string> outputFile;
	  boost::tokenizer<boost::char_separator<char> > tokens(bartCommandLine, sep);
	  for (auto tok: tokens)
	       outputFile.push_back(tok);
	  return outputFile.back();
#endif /* 0 */
     }
}

namespace Gadgetron {

     BartGadget::BartGadget() :
	  BaseClass(),
	  dp{}
     {}

     std::vector<size_t> BartGadget::read_BART_hdr(const std::string& filename)
     {
	  std::vector<size_t> DIMS;
	  
	  std::ifstream infile(filename + ".hdr");
	  if (infile)
	  {
	       std::vector<std::string> tokens;
	       std::string line;

	       while (std::getline(infile, line, '\n'))
	       {
		    tokens.push_back(line);
	       }

	       // Parse the dimensions
	       std::stringstream ss(tokens[1]);
	       std::string items;
	       while (getline(ss, items, ' ')) {
		    DIMS.push_back(std::stoi(items));
	       }
	  }
	  else
	  {
	       GERROR("Failed to header of file: %s\n", filename.c_str());
	  }

	  return DIMS;
     }


     std::pair<std::vector<size_t>, std::vector<std::complex<float>>>
     BartGadget::read_BART_files(const std::string& filename)
     {
	  // Load the header file
	  auto DIMS = read_BART_hdr(filename);

	  // Load the cfl file
	  std::ifstream infile(filename + ".cfl", std::ifstream::binary);
	  if (!infile)
	  {
	       GERROR("Failed to open data of file: %s\n", filename.c_str());
	       return std::make_pair(DIMS, std::vector<std::complex<float>>());
	  }
	  else
	  {
	       // First determine data size
	       infile.seekg(0, infile.end);
	       size_t size(infile.tellg());

	       // Now read data from file
	       infile.seekg(0);
	       std::vector<std::complex<float>> Data(size/sizeof(float)/2.);
	       infile.read(reinterpret_cast<char*>(Data.data()), size);

	       return std::make_pair(DIMS, Data);
	  }
     }

     void BartGadget::replace_default_parameters(std::string & str)
     {
	  std::string::size_type pos = 0u;
	  while ((pos = str.find('$', pos)) != std::string::npos)
	  {
	       auto pos_end = str.find(' ', pos);
	       auto pos_diff = pos_end - pos;
	       std::string tmp = str.substr(pos, pos_diff);
	       tmp.erase(0, 1);
	       if (tmp == std::string("recon_matrix_x"))
		    str.replace(pos, pos_diff, std::to_string(dp.recon_matrix_x));
	       else if (tmp == std::string("recon_matrix_y"))
		    str.replace(pos, pos_diff, std::to_string(dp.recon_matrix_y));
	       else if (tmp == std::string("recon_matrix_z"))
		    str.replace(pos, pos_diff, std::to_string(dp.recon_matrix_z));
	       else if (tmp == std::string("FOV_x"))
		    str.replace(pos, pos_diff, std::to_string(dp.FOV_x));
	       else if (tmp == std::string("FOV_y"))
		    str.replace(pos, pos_diff, std::to_string(dp.FOV_y));
	       else if (tmp == std::string("FOV_z"))
		    str.replace(pos, pos_diff, std::to_string(dp.FOV_z));
	       else if (tmp == std::string("acc_factor_PE1"))
		    str.replace(pos, pos_diff, std::to_string(dp.acc_factor_PE1));
	       else if (tmp == std::string("acc_factor_PE2"))
		    str.replace(pos, pos_diff, std::to_string(dp.acc_factor_PE2));
	       else if (tmp == std::string("reference_lines_PE1"))
		    str.replace(pos, pos_diff, std::to_string(dp.reference_lines_PE1));
	       else if (tmp == std::string("reference_lines_PE2"))
		    str.replace(pos, pos_diff, std::to_string(dp.reference_lines_PE2));
	       else {
		    GERROR( "Unknown default parameter, please see the complete list of available parameters...");
	       }
	       pos = pos_end;
	  }
     }


     bool BartGadget::call_BART(std::string cmdline)
     {
	  GDEBUG_STREAM("Executing BART command: " << cmdline);
	  enum { MAX_ARGS = 256 };

	  // tokenize the command string into argc/argv (basic, hopefully enough)
	  int argc(0);
	  char* argv[MAX_ARGS];

	  auto cmdline_s = std::make_unique<char[]>(cmdline.size()+1);
	  strcpy(cmdline_s.get(), cmdline.c_str());

	  char *p2 = strtok(cmdline_s.get(), " ");
	  while (p2 && argc < MAX_ARGS-1)
	  {
	       argv[argc++] = p2;
	       p2 = strtok(0, " ");
	  }
	  argv[argc] = 0; // FIXME: what is this for ???
	  
	  // boost::char_separator<char> sep(" ");
	  // boost::tokenizer<boost::char_separator<char>> tokens(cmdline.begin(),
	  // 						       cmdline.end(),
	  // 						       sep);
	  // std::vector<std::unique_ptr<char[]>> tmp;
	  // auto k(0UL);
	  // for (auto tok: tokens) {
	  //      tmp.push_back(std::make_unique<char[]>(tok.size()+1));
	  //      strcpy(tmp.back().get(), tok.c_str());
	  //      argv[k++] = tmp.back().get();
	  // }
	  // argc = tmp.size();

	  char out_str[512] = {'\0'};
	  auto ret(in_mem_bart_main(argc, argv, out_str));
	  if (ret == 0) {
	       if (strlen(out_str) > 0) {
		    GINFO(out_str);
	       }
	       return true;
	  }
	  else {
	       GERROR_STREAM("BART command failed with return code: " << ret);
	       return false;
	  }
     }
     
	
     int BartGadget::process_config(ACE_Message_Block * mb)
     {
	  GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

	  /** Let's get some information about the incoming data **/
	  ISMRMRD::IsmrmrdHeader h;
	  try
	  {
	       deserialize(mb->rd_ptr(), h);
	  }
	  catch (...)
	  {
	       GDEBUG("BartGadget::process_config: Failed to parse incoming ISMRMRD Header");
	  }


	  for (const auto& enc: h.encoding)
	  {
	       auto recon_space = enc.reconSpace;

	       GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Encoding matrix size: " << recon_space.matrixSize.x << " " << recon_space.matrixSize.y << " " << recon_space.matrixSize.z);
	       GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Encoding field_of_view : " << recon_space.fieldOfView_mm.x << " " << recon_space.fieldOfView_mm.y << " " << recon_space.fieldOfView_mm.z);
	       dp.recon_matrix_x = recon_space.matrixSize.x;
	       dp.recon_matrix_y = recon_space.matrixSize.y;
	       dp.recon_matrix_z = recon_space.matrixSize.z;
	       dp.FOV_x = static_cast<uint16_t>(recon_space.fieldOfView_mm.x);
	       dp.FOV_y = static_cast<uint16_t>(recon_space.fieldOfView_mm.y);
	       dp.FOV_z = static_cast<uint16_t>(recon_space.fieldOfView_mm.z);

	       if (!enc.parallelImaging)
	       {
		    GDEBUG_STREAM("BartGadget::process_config: Parallel Imaging not enable...");
	       }
	       else
	       {
		    auto p_imaging = *h.encoding.front().parallelImaging;
		    GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: acceleration Factor along PE1 is " << p_imaging.accelerationFactor.kspace_encoding_step_1);
		    GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: acceleration Factor along PE2 is " << p_imaging.accelerationFactor.kspace_encoding_step_2);
		    dp.acc_factor_PE1 = p_imaging.accelerationFactor.kspace_encoding_step_1;
		    dp.acc_factor_PE2 = p_imaging.accelerationFactor.kspace_encoding_step_2;

		    if (p_imaging.accelerationFactor.kspace_encoding_step_2 > 1) {
			 GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Limits of the size of the calibration region (PE1) " << h.userParameters->userParameterLong[0].name << " is " << h.userParameters->userParameterLong[0].value);
			 GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Limits of the size of the calibration region (PE2) " << h.userParameters->userParameterLong[1].name << " is " << h.userParameters->userParameterLong[1].value);
			 dp.reference_lines_PE1 = h.userParameters->userParameterLong[0].value;
			 dp.reference_lines_PE2 = h.userParameters->userParameterLong[1].value;
		    }
		    else if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1) {
			 GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Limits of the size of the calibration region (PE1) " << h.userParameters->userParameterLong[0].name << " is " << h.userParameters->userParameterLong[0].value);
			 dp.reference_lines_PE1 = h.userParameters->userParameterLong[0].value;
		    }

		    auto calib = *p_imaging.calibrationMode;
		    auto separate = (calib.compare("separate") == 0);
		    auto embedded = (calib.compare("embedded") == 0);
		    auto external = (calib.compare("external") == 0);
		    auto interleaved = (calib.compare("interleaved") == 0);
		    auto other = (calib.compare("other") == 0);

		    if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1 || p_imaging.accelerationFactor.kspace_encoding_step_2 > 1)
		    {
			 if (interleaved) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode INTERLEAVE ");
			 }
			 else if (embedded) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode EMBEDDED");
			 }
			 else if (separate) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode SEPERATE");
			 }
			 else if (external) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode EXTERNAL");
			 }
			 else if (other) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode OTHER");
			 }
			 else {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Something went terribly wrong, this should never happen!");
			      return GADGET_FAIL;
			 }
		    }

	       }

	  }
	  return GADGET_OK;
     }

     int BartGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
     {        
	  static std::mutex mtx;
	  std::lock_guard<std::mutex> guard(mtx);

	  // Check status of bart commands script
	  std::string CommandScript = AbsoluteBartCommandScript_path.value() + "/" + BartCommandScript_name.value();
	  if (!boost::filesystem::exists(CommandScript))
	  {
	       GERROR("Can't find bart commands script: %s!\n", CommandScript.c_str());
	       return GADGET_FAIL;
	  }

	  // set the permission for the script
#ifdef _WIN32
	  try
	  {
	       boost::filesystem::permissions(CommandScript, boost::filesystem::all_all);
	  }
	  catch (...)
	  {
	       GERROR("Error changing the permission of the command script.\n");
	       return GADGET_FAIL;
	  }
#else
	  // in case an older version of boost is used in non-win system
	  // the system call is used
	  int res = chmod(CommandScript.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH);
	  if (res != 0)
	  {
	       GERROR("Error changing the permission of the command script.\n");
	       return GADGET_FAIL;
	  }
#endif // _WIN32

	  // Check status of the folder containing the generated files (*.hdr & *.cfl)
	  
	  if (BartWorkingDirectory_path.value().empty()){
	       GERROR("Error: No BART working directory provided!");
	       return GADGET_FAIL;
	  }

	  time_t rawtime;
	  char buff[80];
	  time(&rawtime);
	  strftime(buff, sizeof(buff), "%H_%M_%S__", localtime(&rawtime));
	  std::mt19937::result_type seed = static_cast<unsigned long>(time(0));
	  auto dice_rand = std::bind(std::uniform_int_distribution<int>(1, 10000), std::mt19937(seed));
	  std::string time_id(buff + std::to_string(dice_rand()));
	  // Get the current process ID
	  std::string threadId = boost::lexical_cast<std::string>(boost::this_thread::get_id());
	  unsigned long threadNumber = 0;
	  sscanf(threadId.c_str(), "%lx", &threadNumber);
	  std::string outputFolderPath = BartWorkingDirectory_path.value() + "bart_" + time_id + "_" + std::to_string(threadNumber);
	  std::string generatedFilesFolder = outputFolderPath;
                
	  boost::filesystem::path dir(generatedFilesFolder);
	  if (!boost::filesystem::exists(dir) || !boost::filesystem::is_directory(dir)){
	       if (boost::filesystem::create_directories(dir)){
		    GDEBUG("Folder to store *.hdr & *.cfl files is %s\n", generatedFilesFolder.c_str());
	       }
	       else {
		    GERROR("Folder to store *.hdr & *.cfl files doesn't exist...\n");
		    return GADGET_FAIL;
	       }
	  }

	  generatedFilesFolder += "/";

	  internal::ScopeGuard cleanup_guard(outputFolderPath);

	  /*USE WITH CAUTION*/
	  if (boost::filesystem::exists(generatedFilesFolder) && isBartFolderBeingCachedToVM.value() && !isBartFileBeingStored.value())
	  {
	       std::ostringstream cmd;
	       cmd << "mount -t tmpfs -o size" << AllocateMemorySizeInMegabytes.value() << "M, mode=0755 tmpfs " << generatedFilesFolder;
	       if (system(cmd.str().c_str())) {
		    return GADGET_FAIL;
	       }
	  }

	  std::vector<long> DIMS_ref, DIMS;

	  /*** WRITE REFERENCE AND RAW DATA TO FILES ***/

	  for (auto recon_bit: m1->getObjectPtr()->rbit_)
	  {
	       // Grab a reference to the buffer containing the reference data
	       hoNDArray< std::complex<float> >& data_ref = (*recon_bit.ref_).data_;
	       // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
	       DIMS_ref.emplace_back(data_ref.get_size(0));
	       DIMS_ref.emplace_back(data_ref.get_size(1));
	       DIMS_ref.emplace_back(data_ref.get_size(2));
	       DIMS_ref.emplace_back(data_ref.get_size(3));
	       DIMS_ref.emplace_back(data_ref.get_size(4));
	       DIMS_ref.emplace_back(data_ref.get_size(5));
	       DIMS_ref.emplace_back(data_ref.get_size(6));

	       // Grab a reference to the buffer containing the image data
	       auto& data = recon_bit.data_.data_;
	       // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
	       DIMS.emplace_back(data.get_size(0));
	       DIMS.emplace_back(data.get_size(1));
	       DIMS.emplace_back(data.get_size(2));
	       DIMS.emplace_back(data.get_size(3));
	       DIMS.emplace_back(data.get_size(4));
	       DIMS.emplace_back(data.get_size(5));
	       DIMS.emplace_back(data.get_size(6));

	       /* The reference data will be pointing to the image data if there is
		  no reference scan. Therefore, we won't write the reference data
		  into files if it's pointing to the raw data.*/
	       if (DIMS_ref != DIMS) {
		    // write_BART_Files(std::string(generatedFilesFolder + "meas_gadgetron_ref"), DIMS_ref, data_ref);
		    register_mem_cfl_non_managed("meas_gadgetron_ref", DIMS.size(), &DIMS_ref[0], &data_ref[0]);
	       }

	       // write_BART_Files(std::string(generatedFilesFolder + "meas_gadgetron"), DIMS, data);
	       register_mem_cfl_non_managed("meas_gadgetron", DIMS.size(), &DIMS[0], &data[0]);
	  }

	  /* Before calling Bart let's do some bookkeeping */
	  std::ostringstream cmd1, cmd2, cmd3;
	  std::replace(generatedFilesFolder.begin(), generatedFilesFolder.end(), '\\', '/');

	  if (DIMS_ref != DIMS)
	  {
	       cmd1 << "bart resize -c 0 " << DIMS[0] << " 1 " << DIMS[1] << " 2 " << DIMS[2] << " meas_gadgetron_ref reference_data";
	       // Pass commands to Bart
	       if (!call_BART(cmd1.str())) {
		    return GADGET_FAIL;
	       }
	  }

	  if (DIMS[4] != 1)
	       cmd2 << "bart reshape 1023 " << DIMS[0] << " " << DIMS[1] << " " << DIMS[2] << " " << DIMS[3] << " 1 1 1 " << DIMS[5] << " " << DIMS[6] << " " << DIMS[4] << " meas_gadgetron input_data";
	  else	
	       cmd2 << "bart scale 1.0 meas_gadgetron input_data";

	  if (!call_BART(cmd2.str())) {
	       return GADGET_FAIL;
	  }

	  /*** CALL BART COMMAND LINE from the scripting file***/

	  std::string Line, Commands_Line;
	  std::fstream inputFile(CommandScript);
	  if (inputFile.is_open())
	  {
	       while (getline(inputFile, Line))
	       {
		    // crop comment
		    Line = Line.substr(0, Line.find_first_of("#"));

		    internal::trim(Line);
		    if (Line.empty() || Line.compare(0, 4, "bart") != 0)
			 continue;
				
		    replace_default_parameters(Line);
		    GDEBUG("%s\n", Line.c_str());

		    if (!call_BART(Line)) {
			 return GADGET_FAIL;
		    }

		    // system(std::string("cd " + generatedFilesFolder + "&&" + Line).c_str());
				
		    Commands_Line = Line;
	       }
	       inputFile.close();
	  }
	  else
	  {
	       GERROR("Unable to open %s\n", CommandScript.c_str());
	       return GADGET_FAIL;
	  }

	  std::string outputFile = internal::get_output_filename(Commands_Line);
	  std::string outputFileReshape = outputFile + "_reshape";

	  // Reformat the data back to gadgetron format
	  std::vector<long> header(16);
	  load_mem_cfl(outputFile.c_str(), header.size(), header.data());
	  // auto header = read_BART_hdr(generatedFilesFolder + outputFile);
	  cmd3 << "bart reshape 1023 " << header[0] << " " << header[1] << " " << header[2] << " " << header[3] << " " << header[9] * header[4] <<
	       " " << header[5] << " " << header[6] << " " << header[7] << " " << header[8] << " 1 " << outputFile << " " << outputFileReshape;

	  if (!call_BART(cmd3.str())) {
	       return GADGET_FAIL;
	  }
	  
	  /**** READ FROM BART FILES ***/
	  std::vector<long> DIMS_OUT(16);
	  std::complex<float>* data(reinterpret_cast<std::complex<float>*>(load_mem_cfl(outputFileReshape.c_str(), DIMS_OUT.size(), DIMS_OUT.data())));

	  if (data == 0 || data == nullptr) {
	       GERROR("Failed to retrieve data from in-memory CFL file!");
	       return GADGET_FAIL;
	  }

	  // std::vector<std::size_t> DIMS_OUT;
	  // std::vector<std::complex<float>> data;
	  // std::tie(DIMS_OUT, data) = read_BART_files(generatedFilesFolder + outputfileReshape);

	  if (!isBartFileBeingStored.value())
	       internal::cleanup(outputFolderPath);
	  else
	       cleanup_guard.dismiss();

	  auto ims = std::make_unique<GadgetContainerMessage<IsmrmrdImageArray>>();
	  IsmrmrdImageArray & imarray = *ims->getObjectPtr();

	  // Grab data from BART files
	  std::vector<size_t> BART_DATA_dims(1);
	  BART_DATA_dims[0] = std::accumulate(DIMS_OUT.begin(), DIMS_OUT.end(), 1, std::multiplies<size_t>());
	  hoNDArray<std::complex<float>> DATA(BART_DATA_dims, data);
	  // hoNDArray<std::complex<float>> DATA(BART_DATA_dims, &data[0]);

	  // The image array data will be [E0,E1,E2,1,N,S,LOC]
	  std::vector<size_t> data_dims(7);
	  std::copy(DIMS_OUT.begin(), DIMS_OUT.begin()+7, data_dims.begin());

	  DATA.reshape(data_dims);

	  // Extract the first image from each time frame (depending on the number of maps generated by the user)
	  std::vector<size_t> data_dims_Final(7);
	  data_dims_Final[0] = DIMS_OUT[0];
	  data_dims_Final[1] = DIMS_OUT[1];
	  data_dims_Final[2] = DIMS_OUT[2];
	  data_dims_Final[3] = DIMS_OUT[3];
	  assert(header[4] > 0);
	  data_dims_Final[4] = DIMS_OUT[4] / header[4];
	  data_dims_Final[5] = DIMS_OUT[5];
	  data_dims_Final[6] = DIMS_OUT[6];
	  imarray.data_.create(&data_dims_Final);

	  std::vector<std::complex<float> > DATA_Final;
	  DATA_Final.reserve(std::accumulate(data_dims_Final.begin(), data_dims_Final.end(), 1, std::multiplies<size_t>()));

	  for (uint16_t loc = 0; loc < data_dims[6]; ++loc) {
	       for (uint16_t s = 0; s < data_dims[5]; ++s) {
		    for (uint16_t n = 0; n < data_dims[4]; n += header[4]) {

			 //Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
			 //Each chunk will be [E0,E1,E2,CHA] big
			 std::vector<size_t> chunk_dims(4), Temp_one_1d(1);
			 chunk_dims[0] = data_dims_Final[0];
			 chunk_dims[1] = data_dims_Final[1];
			 chunk_dims[2] = data_dims_Final[2];
			 chunk_dims[3] = data_dims_Final[3];
			 hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &DATA(0, 0, 0, 0, n, s, loc));

			 Temp_one_1d[0] = chunk_dims[0] * chunk_dims[1] * chunk_dims[2] * chunk_dims[3];
			 chunk.reshape(Temp_one_1d);
			 DATA_Final.insert(DATA_Final.end(), chunk.begin(), chunk.end());
		    }
	       }
	  }

	  std::copy(DATA_Final.begin(), DATA_Final.end(), imarray.data_.begin());

	  // Fill image header 
	  for (size_t it = 0, rbit_size = m1->getObjectPtr()->rbit_.size(); it != rbit_size; ++it)
	  {
	       compute_image_header(m1->getObjectPtr()->rbit_[it], imarray, it);
	       send_out_image_array(m1->getObjectPtr()->rbit_[it], imarray, it, image_series.value() + (static_cast<int>(it) + 1), GADGETRON_IMAGE_REGULAR);
	  }

	  m1->release();
	  return GADGET_OK;
     }

     GADGET_FACTORY_DECLARE(BartGadget)
}
