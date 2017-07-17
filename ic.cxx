#include <iostream>
#include <fstream>      
#include <string> 
#include <ctime>    // added by khurram
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "linterp.h"
#include "ic-read-temp.hpp"
#include "ic.hpp"
using namespace std;    //added by khurram



namespace {

    class Zone
    {
    public:
        virtual ~Zone() {};
        virtual bool contains(const double x[NDIMS]) const = 0;
    };

    class Empty_zone : public Zone
    {
    public:
        bool contains(const double x[NDIMS]) const {return false;}
    };


    class Planar_zone : public Zone
    {
    private:
        const double az, incl;
        const double halfwidth; // in meter
#ifdef THREED
        const double ymin, ymax; // in meter
#endif
        const double zmin, zmax; // in meter
        const double *x0;

    public:
        Planar_zone(const double center[NDIMS], double azimuth, double inclination, double halfwidth_,
#ifdef THREED
                    double ymin_, double ymax_,
#endif
                    double zmin_, double zmax_) :
            az(std::tan(azimuth * DEG2RAD)), incl(1/std::tan(inclination * DEG2RAD)), halfwidth(halfwidth_),
#ifdef THREED
            ymin(ymin_), ymax(ymax_),
#endif
            zmin(zmin_), zmax(zmax_),
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {}

        bool contains(const double x[NDIMS]) const
        {
            // Is x within halfwidth distance to a plane cutting through x0?
            return (x[NDIMS-1] > zmin &&
                    x[NDIMS-1] < zmax &&
#ifdef THREED
                    x[1] > ymin &&
                    x[1] < ymax &&
#endif
                    std::fabs( (x[0] - x0[0])
#ifdef THREED
                               - az * (x[1] - x0[1])
#endif
                               + incl * (x[NDIMS-1] - x0[NDIMS-1]) ) < halfwidth );
        }
    };


    class Ellipsoidal_zone : public Zone
    {
    private:
        const double *x0;
        double semi_axis2[NDIMS];

    public:
        Ellipsoidal_zone(const double center[NDIMS], const double semi_axis[NDIMS]) :
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {
            for(int i=0; i<NDIMS; i++)
                semi_axis2[i] =  semi_axis[i] * semi_axis[i];
        }

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])/semi_axis2[0]
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])/semi_axis2[1]
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])/semi_axis2[NDIMS-1] < 1 );
        }
    };

} // anonymous namespace


void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, tensor_t &strain,
                          double &compensation_pressure)
// void initial_stress_state(const Param &param, const Variables &var,
//                           tensor_t &stress, double_vec &stressyy, tensor_t &strain,
//                           double &compensation_pressure, double_vec &stressxx, double_vec &stressxz, double_vec &stresszz)

{
    // if (param.control.gravity == 0) {
    //     compensation_pressure = 0;
    //     return;
    // }



    //------------------------------------------
    //loading the input x and y vectors 
    //-----------------------------------------
    int rows= 1601;    //y grid points
    int cols=801;
    vector<double> x_vector(rows);
    vector<double> y_vector(cols);

    ifstream x_coord("x_vector.in");
    if (!x_coord) 
    {
        cout << "Cannot open x_vector file.\n";
        abort();
    }
    
    for (int i= 0; i < rows; i++)
    {
        x_coord >> x_vector[i];
        x_vector[i]= x_vector[i] *1000.0; //changing values from km to m
    } 
    x_coord.close(); 



    ifstream y_coord("y_vector.in");
    if (!y_coord) 
    {
        cout << "Cannot open y_vector file.\n";
        abort();
    }
    
    for (int i= 0; i < cols; i++)
    {
        y_coord >> y_vector[i];
        y_vector[i]= y_vector[i] *1000.0; //changing values from km to m
    }
    y_coord.close(); 
    //----------------------------------------------
    
    // Loading the input stress 
    //------------------------------------------

    // initialize stresses
    double** stress_sxx = new double*[rows];
    double** stress_sxy = new double*[rows];
    double** stress_syy = new double*[rows];
    for (int i = 0; i < rows; ++i)
    {
        stress_sxx[i] = new double[cols];
        stress_sxy[i] = new double[cols];
        stress_syy[i] = new double[cols];
    }

    
    // loading sxx stress component
    ifstream stress_sxx_file("sxx_stress.in");
    if (!stress_sxx_file) 
    {
        cout << "Cannot open sxx_stress file.\n";
        abort();
    }
    
    for (int i = 0; i < rows; i++) 
    {
        for (int j = 0; j < cols; j++) 
        {
            stress_sxx_file >> stress_sxx[i][j];
            stress_sxx[i][j]= stress_sxx[i][j] * 1000000.0 ; // convert it into pascal    
        }
    }
    stress_sxx_file.close();

        // loading sxy stress component
    ifstream stress_sxy_file("sxy_stress.in");
    if (!stress_sxy_file) 
    {
        cout << "Cannot open sxy_stress file.\n";
        abort();
    }
    
    for (int i = 0; i < rows; i++) 
    {
        for (int j = 0; j < cols; j++) 
        {
            stress_sxy_file >> stress_sxy[i][j];
            stress_sxy[i][j]= stress_sxy[i][j] * 1000000.0;     
        }
    }
    stress_sxy_file.close();

        // loading syy stress component
    ifstream stress_syy_file("syy_stress.in");
    if (!stress_syy_file) 
    {
        cout << "Cannot open syy_stress file.\n";
        abort();
    }
    
    for (int i = 0; i < rows; i++) 
    {
        for (int j = 0; j < cols; j++) 
        {
            stress_syy_file >> stress_syy[i][j];
            stress_syy[i][j] = stress_syy[i][j] * 1000000.0 ;    
        }
    }
    stress_syy_file.close();


    //--------------------------------------------------
    // Interpolation part 
    //--------------------------------------------------

    // note that we will pass in a sequence of iterators pointing to the beginning of each grid
    std::vector< std::vector<double>::iterator > grid_iter_list;
    grid_iter_list.push_back(x_vector.begin());
    grid_iter_list.push_back(y_vector.begin());
    // the size of the grid in each dimension
    array<int,2> grid_sizes;
    grid_sizes[0] = rows;
    grid_sizes[1] = cols;
  
  // total number of elements
    int num_elements = grid_sizes[0] * grid_sizes[1];
    

  
  // fill in the values of f(x) at the gridpoints. 
  // we will pass in a contiguous sequence, values are assumed to be laid out C-style
   std::vector<double> f_values_sxx(num_elements);
   std::vector<double> f_values_sxy(num_elements);
   std::vector<double> f_values_syy(num_elements);

    for (int i=0; i<grid_sizes[0]; i++) {
         for (int j=0; j<grid_sizes[1]; j++) {
            f_values_sxx[i*grid_sizes[1] + j] = stress_sxx[i][j];
            f_values_sxy[i*grid_sizes[1] + j] = stress_sxy[i][j];
            f_values_syy[i*grid_sizes[1] + j] = stress_syy[i][j];
        }
    }




    // construct the interpolator. the last two arguments are pointers to the underlying data
    
    InterpMultilinear<2, double> interp_ML_sxx(grid_iter_list.begin(), grid_sizes.begin(), f_values_sxx.data(), f_values_sxx.data() + num_elements);
    InterpMultilinear<2, double> interp_ML_sxy(grid_iter_list.begin(), grid_sizes.begin(), f_values_sxy.data(), f_values_sxy.data() + num_elements);
    InterpMultilinear<2, double> interp_ML_syy(grid_iter_list.begin(), grid_sizes.begin(), f_values_syy.data(), f_values_syy.data() + num_elements);
    
    // interpolate one value
    
    array<double,2> args = {45000.0, 22000.0};
    printf("%f, %f stressyy-> %f\n", args[0], args[1], interp_ML_syy.interp(args.begin()));



    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    // --------------------------------------------------------------------
    // Assign the values of stresses to the berycenter in each element
    // --------------------------------------------------------------------

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        double xcenter = 0;

        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
            xcenter+=  (*var.coord)[conn[i]][0];   //this line is added by khurram
        }
        zcenter /= NODES_PER_ELEM;
        xcenter /= NODES_PER_ELEM;
        array<double,2> args = {xcenter, zcenter};
        //stressyy[e]= interp_ML_syy.interp(args.begin());   //need to change it to stresszz
        // stressxz[e]= interp_ML_sxy.interp(args.begin()); 
        // stressxx[e]= interp_ML_sxx.interp(args.begin());


        double p = ref_pressure(param, zcenter);
        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
        }

        stress[e][0] = interp_ML_sxx.interp(args.begin());
        stress[e][1] = interp_ML_syy.interp(args.begin());
        stress[e][2] = interp_ML_sxy.interp(args.begin());
        
        for (int i=0; i<NDIMS; ++i) {
            //stress[e][i] = -p;  // in 2D. i=0 means 'xx'; i=1 means 'zz'; i=2 (==NDIMS) means 'xz'. 
            

            strain[e][i] = -p / ks / NDIMS;
        }
    //     if (param.mat.is_plane_strain)
    //         stressyy[e] = -p;
    }

    //compensation_pressure = ref_pressure(param, -param.mesh.zlength);

    // freeing the matrices allocated at the top of the function
    for (int i = 0; i < rows; ++i)
    {
        delete [] stress_sxx[i];
        delete [] stress_sxy[i];
        delete [] stress_syy[i];
    }    
    delete [] stress_sxx;
    delete [] stress_sxy;
    delete [] stress_syy;

}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    Zone *weakzone;

    // TODO: adding different types of weak zone
    double plane_center[NDIMS]; // this variable must outlive weakzone
    switch (param.ic.weakzone_option) {
    case 0:
        weakzone = new Empty_zone();
        break;
    case 1:
        // a planar weak zone, cut through top center
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Planar_zone(plane_center,
                                   param.ic.weakzone_azimuth,
                                   param.ic.weakzone_inclination,
                                   param.ic.weakzone_halfwidth * param.mesh.resolution,
#ifdef THREED
                                   param.ic.weakzone_y_min * param.mesh.ylength,
                                   param.ic.weakzone_y_max * param.mesh.ylength,
#endif
                                   -param.ic.weakzone_depth_max * param.mesh.zlength,
                                   -param.ic.weakzone_depth_min * param.mesh.zlength);
        break;
    case 2:
        // a ellipsoidal weak zone
        double semi_axis[NDIMS];
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
        semi_axis[0] = param.ic.weakzone_xsemi_axis;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
        semi_axis[1] = param.ic.weakzone_ysemi_axis;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        semi_axis[NDIMS-1] = param.ic.weakzone_zsemi_axis;
        weakzone = new Ellipsoidal_zone(plane_center, semi_axis);
        break;
    default:
        std::cerr << "Error: unknown weakzone_option: " << param.ic.weakzone_option << '\n';
        std::exit(1);
    }

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[NDIMS] = {0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }

        if (weakzone->contains(center))
            plstrain[e] = param.ic.weakzone_plstrain;

        // Find the most abundant marker mattype in this element
        // int_vec &a = (*var.elemmarkers)[e];
        // int material = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    }

    delete weakzone;
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    switch(param.ic.temperature_option) {
    case 0:
        {
            const double age = param.ic.oceanic_plate_age_in_yr * YEAR2SEC;
            const MatProps &mat = *var.mat;
            const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element

            for (int i=0; i<var.nnode; ++i) {
                double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * age);
                temperature[i] = param.bc.surface_temperature +
                    (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
            }
            break;
        }
    case 90:
        read_external_temperature_from_comsol(param, var, *var.temperature);
        break;
    default:
        std::cout << "Error: unknown ic.temperature option: " << param.ic.temperature_option << '\n';
        std::exit(1);
    }
}


