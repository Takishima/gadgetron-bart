/****************************************************************************************************************************
 * Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Author: Mahamadou Diakite, PhD.
 * Institution: National Institutes of Health (NIH)
 * Lang: C++
 * Date: 10/15/2017
 * Version: 1.0.0
 ****************************************************************************************************************************/

#ifndef BART_GADGET_H
#define BART_GADGET_H

#include "GenericReconGadget.h"
#include "gadgetron_mricore_export.h"
#include "mri_core_data.h"
#include <boost/filesystem.hpp>
#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>
#include "gadgetron_home.h"

#if defined (WIN32)
#ifdef __BUILD_GADGETRON_bartgadget__
#define EXPORTGADGETS_bartgadget __declspec(dllexport)
#else
#define EXPORTGADGETS_bartgadget __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_bartgadget
#endif


namespace Gadgetron {
	
     // The user is free to add more parameters as the need arises.
     struct Default_parameters 
     {
	  uint16_t recon_matrix_x;
	  uint16_t recon_matrix_y;
	  uint16_t recon_matrix_z;
	  uint16_t FOV_x;
	  uint16_t FOV_y;
	  uint16_t FOV_z;
	  uint16_t acc_factor_PE1;
	  uint16_t acc_factor_PE2;
	  uint16_t reference_lines_PE1;
	  uint16_t reference_lines_PE2;
     };

     class EXPORTGADGETS_bartgadget BartGadget final : public GenericReconGadget
     {

     public:
	  GADGET_DECLARE(BartGadget)
		
	       using BaseClass = GenericReconGadget;

	  BartGadget();
	  ~BartGadget() = default;

     protected:
	  GADGET_PROPERTY(isVerboseON, bool, "Display some information about the incoming data", false);
	  GADGET_PROPERTY(BartWorkingDirectory_path, std::string, "Absolute path to temporary file location", "/tmp/gadgetron/");
	  GADGET_PROPERTY(AbsoluteBartCommandScript_path, std::string, "Absolute path to bart script(s)", get_gadgetron_home().string() + "/share/gadgetron/bart");
	  GADGET_PROPERTY(BartCommandScript_name, std::string, "Script file containing BART command(s) to be loaded", "");
	  GADGET_PROPERTY(isBartFileBeingStored, bool, "Store BART file on the disk", false);
	  GADGET_PROPERTY(image_series, int, "Set image series", 0);

	  /*Caution: this option must be enable only if the user has root privilege and able to allocation virtual memory*/
	  GADGET_PROPERTY(isBartFolderBeingCachedToVM, bool, "Mount bart directory to virtual memory (tmpfs) for better performance", false);
	  GADGET_PROPERTY(AllocateMemorySizeInMegabytes, int, "Allocate memory to bart directory", 50);

	  int process_config(ACE_Message_Block* mb);
	  int process(GadgetContainerMessage<IsmrmrdReconData>* m1);		

     private:
	  Default_parameters dp;
		
	  void replace_default_parameters(std::string &str);
	  
	  bool call_BART(std::string cmdline);	  
     };

     // Read BART files
     std::vector<size_t> read_BART_hdr(const std::string& filename);
     std::pair< std::vector<size_t>, std::vector<std::complex<float> > > read_BART_files(const std::string& filename);


     template<typename int_t>
     void write_BART_hdr(std::string filename, const std::vector<int_t>& DIMS)
     {
	  constexpr size_t MAX_DIMS = 16;
	  std::vector<size_t> v(MAX_DIMS, 1);
	  assert(DIMS.size() < MAX_DIMS);
	  std::copy(DIMS.cbegin(), DIMS.cend(), v.begin());
	  std::ofstream pFile(filename + ".hdr");
	  if (!pFile)
	       GERROR("Failed to write into header file: %s\n", filename);
	  pFile << "# Dimensions\n";
	  std::copy(v.cbegin(), v.cend(), std::ostream_iterator<size_t>(pFile, " "));
     }

     template<typename int_t>
     void write_BART_Files(std::string filename, const std::vector<int_t>& DIMS, const std::vector<std::complex<float>>& DATA)
     {
	  write_BART_hdr(filename, DIMS);
	  std::ofstream pFile(filename + ".cfl", std::ofstream::out | std::ofstream::binary);
	  if (!pFile)
	       GERROR("Failed to write into CFL file: %s\n", filename);
	  pFile.write(reinterpret_cast<const char*>(&DATA[0]), DATA.size() * sizeof(std::complex<float>));
     }
     
     template<typename int_t>
     void write_BART_Files(std::string filename, const std::vector<int_t>&DIMS, const hoNDArray<std::complex<float>>& DATA)
     {
	  write_BART_hdr(filename, DIMS);
	  std::ofstream pFile(filename + ".cfl", std::ofstream::out | std::ofstream::binary);
	  if (!pFile)
	       GERROR("Failed to write into CFL file: %s\n", filename);
	  pFile.write(reinterpret_cast<const char*>(&DATA[0]), DATA.get_number_of_bytes());
     }
} // namespace Gadgetron
#endif //BART_GADGET_H
