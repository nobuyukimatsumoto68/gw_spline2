#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include <Eigen/Dense>




#include "header.h"





int main(int argc, char* argv[]){
  // int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------

  std::string mass = "0p4000";
  if(argc>=2) mass = argv[1];

  int ibetamin, ibetamax;
  if(argc>=3) ibetamin = atoi(argv[2]);
  if(argc>=4) ibetamax = atoi(argv[3]);

  const int Nt=8;
  int Ns=32;
  if(argc>=5) Ns = atoi(argv[4]);

  int nptsx = 40;
  if(argc>=6) nptsx = atoi(argv[5]);
  int nptsy = nptsx;



  const std::string basedir = "/mnt/hdd_barracuda/llnl/reweight/data/"+std::to_string(Ns)+"b_v14/";


  const double minOx = -0.2;
  const double maxOx = 0.6;

  const int irow1 = nptsy/2-1;



  std::vector<double> xpts(nptsx);
  const double deltax = (maxOx-minOx)/nptsx;
  for(int iptx=0; iptx<nptsx; iptx++) xpts[iptx] = minOx + deltax*(iptx+0.5);


  std::vector<double> betas;
  {
    const std::string filename = basedir + "/m"+mass+"avghist_ibx"+std::to_string(0)+"_iby"+std::to_string(0)+"_nojkmeas.bin";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    betas = f.getDataSet("beta").read<std::vector<double>>();
  }
  const int nbeta_meas = betas.size();


  Eigen::ArrayXXd Gamma_list;
  get_Gamma( Gamma_list, nbeta_meas, nptsx, nptsy, irow1, basedir, mass, Ns, Nt );


  for(int ibeta=ibetamin; ibeta<=ibetamax; ibeta++){
    std::cout << "# ibeta: " << ibeta << std::endl;
    // Eigen::ArrayXd yy = Gamma_list.matrix().block( ibeta, 0, 1,nptsx ).array();
    // const Eigen::MatrixXd tmp = Gamma_list.matrix();
    // Eigen::VectorXd yy = tmp.block( ibeta, 0, 1, nptsx );
    // const Eigen::Matrix<double, 1,nptsx> yy = Gamma_list.block( ibeta, 0, 1, nptsx );
    Eigen::MatrixXd yy(1, nptsx);
    yy = Gamma_list.block( ibeta, 0, 1, nptsx );
    // const Eigen::ArrayXd yy = Gamma_list.block<1,nptsx>( ibeta, 0 );
    // for(int i=0; i<xpts.size(); i++){
    //   std::cout << xpts[i] << " " << yy(0,i) << std::endl;
    // }

    // interpolate
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nptsx);
    gsl_spline_init(spline, xpts.data(), yy.data(), nptsx);

    auto Veff = [&](const double x) { return gsl_spline_eval(spline, x, acc); };
    auto Veff_prime = [&](const double x) { return gsl_spline_eval_deriv(spline, x, acc); };
    auto Veff_prime2 = [&](const double x) { return gsl_spline_eval_deriv2(spline, x, acc); };

    // std::cout << "debug. 1" << std::endl;

    double xA, xB;
    get_local_minima( xA, xB, xpts, deltax, Veff_prime, Veff_prime2);
    auto Veff_sub = [&](const double x) { return Veff(x) - Veff(xB); };

    // std::cout << "debug. 2" << std::endl;

    {
      const std::string dir = "./fit_params_"+std::to_string(Ns)+"c_m"+mass+"/";
      std::filesystem::create_directory( dir );

      {
        const int xsize = 1000;
        const double xmin= *std::min_element(xpts.begin(), xpts.end()); // *std::min(xpts); // -0.1;
        const double xmax= *std::max_element(xpts.begin(), xpts.end()); // *std::min(xpts); // -0.1;
        std::vector<double> xs(xsize+1);
        for(int i=0; i<=xsize; i++) xs[i] = (xmax-xmin)/xsize * i + xmin;
        // std::vector<double> xs(xsize-1);
        // for(int i=0; i<xsize-1; i++) xs[i] = (xmax-xmin)/xsize * (i+1) + xmin;
        std::vector<double> Vs;
        for(const double x : xs) Vs.push_back( Veff_sub(x) );

        const std::string out = dir+"Veff_"+std::to_string(ibeta)+".dat";
        write( out, xs, Vs, "VB subtracted" );
      }
      {
        std::vector<double> minima_data{ xA, xB, Veff(xA), Veff(xB) };
        const std::string out = dir+"minima_data_"+std::to_string(ibeta)+".dat";

        std::ofstream ofs(out);
        ofs << std::scientific << std::setprecision(25);
        for(const auto elem : minima_data) ofs << elem << " ";
        ofs << std::endl;
      }
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }


  return 0;
}
