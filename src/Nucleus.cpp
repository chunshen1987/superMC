#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iomanip>

#include "MathBasics.h"
#include "Nucleus.h"
#include "ParameterReader.h"
#include "GaussianNucleonsCal.h"
#include "Box2D.h"
#include "GluonField.h"

using namespace std;

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// --- initializes parameters for WS/Hulthen density distribution of a nucleus
//     function getRandomWS(x,y,z) then returns sampled nucleon coordinate
Nucleus::Nucleus(int a, ParameterReader* paraRdr, int deformed_in, int id)
{
  nuclID = id;

  flag_NN_correlation = paraRdr->getVal("include_NN_correlation");
  GaussianNucleonsCal* gaussCal = new GaussianNucleonsCal(paraRdr);
  double gaussian_entropy_width = gaussCal->width;
  double quark_width = paraRdr->getVal("quark_width");
  // B = sigma_g^2 + (2/3)sigma_q^2. This is when the third quark is determined by CM of nucleon.
  double quark_dist_width = sqrt((3.0/2.0)*(gaussian_entropy_width*gaussian_entropy_width-
                            quark_width*quark_width));

  GFF = (paraRdr->getVal("gluon_field_fluctuations") == 1);
  GluonField::distanceScalingFactor = paraRdr->getVal("gluon_distance_scaling");

  ifstream fin("tables/QuarkPos.txt");
  double x,y,z;
  fin >> x >> y >> z;
  while(!fin.eof()){
    vector<double> temp;
    temp.push_back(x);
    temp.push_back(y);
    temp.push_back(z);
    Particle::quark_pos.push_back(temp);
    fin >> x >> y >> z;
  }
  fin.close();

  Particle::width = gaussian_entropy_width;
  Particle::R = quark_dist_width;
  Quark::width = quark_width;
  delete gaussCal;
  ecm = paraRdr->getVal("ecm");

  // save atomic # of nucleus for later; can be recalled using: getAtomic()
  deformed = deformed_in;
  A = (double)a;
  atomic = a;

  // if nucleon, nothing more needed
  if (a==1) return;

  // generic (crude) default parameterization of radius and surface thickness
  rad=1.12*pow(A,0.333333)-0.86/pow(A,0.333333);
  dr = 0.54;
  density0 = 0.17;

  //taken from C.W.De Jager et al. Atom.Data.Nucl.Data Tabl.36, 495 (1987).
  //if(a==197){
  //  rad = 6.38;
  //  dr = 0.535;
  //  density0 = 0.1695;
  //}else if(a==63){
  //  rad = 4.20641;
  //  dr = 0.5877;
  //  density0 = 0.1686;
  //}

  // to reproduce ordinary Woods-Saxon for finite size nucleon (T. Hirano)
  if(a == 197){
    //Reparametrization (Au)
    rad=6.42;
    dr=0.45;
    density0 = 0.1695;
  }else if(a==63){
    //Reparametrization (Cu)
    rad=4.28;
    dr=0.5;
    density0 = 0.1686;
  }else if(a==238){
    //Reparametrization (U)
    rad=6.86;
    dr=0.44;
    density0 = 0.166;

    //Taken from P.Filip et al.,PRC80,054903(2009)
    //      rad=6.81;
    //      dr=0.54;
    //      density0 = 0.17;
  }else if(a==208){
    //Reparametrization (Pb)
    rad = 6.67;
    dr = 0.44;
    density0 = 0.161;
  }

  // default deformation parameters, 0: no deform, 1: deform (05/03/2010 by TH)
  beta2 = 0.0;
  beta4 = 0.0;
  ctr = 1.0;
  phir = 0.0;
  rmaxCut = rad + 2.5;
  rwMax = 1.0/(1.0+exp(-rad/dr));


  if(deformed){
    // Parameters taken from Moller et al., Atomic Data and Nuclear Data
    // Tables, 59, 185-381, (1995).
    if(atomic == 197){
      //(Au)
      beta2 = -0.13;
      beta4 = -0.03;
    }else if(atomic == 63){
      //(Cu)
      beta2 = 0.162;
      beta4 = 0.006;
    }else if(atomic == 238){
      //(U)
      beta2 = 0.28;
      beta4 = 0.093;
    }else{
      cerr << "Mass number not available for reparametrization" << endl;
    }
  }
  else
  {
    beta2 = 0.0; beta4 = 0.0;
  }

  if (atomic == 3) // read in triton position
  {
    readin_triton_position();
  }

  if(flag_NN_correlation == 1 && (atomic == 197 || atomic == 208))
  {
     readin_nucleon_positions();
  }
}

Nucleus::~Nucleus()
{
   if(flag_NN_correlation == 1 && (atomic == 197 || atomic == 208))
   {
      for(int iconf = 0; iconf < n_configuration; iconf++)
      {
         for(int ia = 0; ia < atomic; ia++)
            delete [] nucleon_pos_array[iconf][ia];
         delete [] nucleon_pos_array[iconf];
      }
      delete [] nucleon_pos_array;
   }
   clearNucleons();
}

void Nucleus::populate(double xCenter, double yCenter) {
    nPart = 0;

    double rmin=0.9*0.9; // minimal nucleon separation (squared; in fm^2).
    double cx, phi;
    double xcm=0.0, ycm=0.0, zcm=0.0;
    cx = 1.0-2.0*drand48();
    phi = 2*M_PI*drand48();
    lastCx = cx;
    lastPh = phi;
    setRotation(cx, phi);
    if (atomic == 1) {
        nucleons.push_back(new Particle(xCenter, yCenter, 0.0)); // shift along x-axis
    } else if (atomic == 2) {
        // in the case of a deuteron (added by Brian Baker)
        double x1, y1, z1, x2, y2, z2;

        GetDeuteronPosition(x1,y1,z1,x2,y2,z2);
        // shift along x-axis
        nucleons.push_back(new Particle(x1+(xCenter),y1+yCenter,z1));
        nucleons.push_back(new Particle(x2+(xCenter),y2+yCenter,z2));
    } else if (atomic == 3) {
        // in the case of 3He
        double x1, x2, x3, y1, y2, y3, z1, z2, z3;
        // Triton data points from Carlson have center of mass (0,0,0).
        GetTritonPosition(x1,y1,z1,x2,y2,z2,x3,y3,z3);
        nucleons.push_back(new Particle(x1+(xCenter),y1+yCenter,z1));
        nucleons.push_back(new Particle(x2+(xCenter),y2+yCenter,z2));
        nucleons.push_back(new Particle(x3+(xCenter),y3+yCenter,z3));
    } else {
        if (flag_NN_correlation == 1 && (atomic == 197 || atomic == 208)) {
            // load nucleon positions from file
            double **nucleon_positions = new double* [atomic];
            for (int inucleon = 0; inucleon < atomic; inucleon++) {
                nucleon_positions[inucleon] = new double [3];
            }
            get_nucleon_position_with_NN_correlation(nucleon_positions);

            double xcm = 0.0, ycm = 0.0, zcm = 0.0;
            for (int ia = 0; ia < atomic; ia++) {
                double x_local = nucleon_positions[ia][0];
                double y_local = nucleon_positions[ia][1];
                double z_local = nucleon_positions[ia][2];
                xcm += x_local;
                ycm += y_local;
                zcm += z_local;
                nucleons.push_back(new Particle(x_local, y_local, z_local));
            }

            for (int ia = 0; ia < atomic; ia++) {
                double x_local = nucleons[ia]->getX() - xcm/atomic + xCenter;
                double y_local = nucleons[ia]->getY() - ycm/atomic + yCenter;
                double z_local = nucleons[ia]->getZ() - zcm/atomic;
                nucleons[ia]->setX(x_local);
                nucleons[ia]->setY(y_local);
                nucleons[ia]->setZ(z_local);
            }

            // clean up
            for(int inucleon = 0; inucleon < atomic; inucleon++) {
                delete [] nucleon_positions[inucleon];
            }
            delete [] nucleon_positions;
        } else {
            // for large nuclei
            // generate a cloud of nucleons around 0,0,0
            // Then move them to the real center
            for(int ia = 0; ia < atomic; ia++) {
                double x,y,z;
                int icon=0;
                int oversamplePosition = 0;
                do {
                    getDeformRandomWS(x, y, z);
                    icon = 0;
                    oversamplePosition = atomic*((int)nucleons.size()/atomic);
                    for(int i=oversamplePosition; i<(int)nucleons.size();i++) {
                        double x1=nucleons[i]->getX();
                        double y1=nucleons[i]->getY();
                        double z1=nucleons[i]->getZ();
                        double r2 = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1);
                        if(r2 < rmin) {
                            icon=1;
                            break;
                        }
                    }
                } while(icon==1);

                xcm +=x;
                ycm +=y;
                zcm +=z;
                nucleons.push_back(new Particle(x,y,z));
            }
            for (int ia=0;ia<atomic;ia++) {
                // shift center of nucleus
                double x = nucleons[ia]->getX() - xcm/atomic + xCenter;
                double y = nucleons[ia]->getY() - ycm/atomic + yCenter;
                double z = nucleons[ia]->getZ() - zcm/atomic;
                nucleons[ia]->setX(x);
                nucleons[ia]->setY(y);
                nucleons[ia]->setZ(z);
            }
        }
    }

    // Sort them by their xLeft bound. This is to optimize collision detection.
    std::sort(nucleons.begin(), nucleons.end(), sortByXLeft);
    if (GFF)
      setGluonFields();
}

void Nucleus::setGluonFields() {
    for(unsigned int i = 0; i < nucleons.size(); i++) {
        nucleons[i]->setGluonField(new GluonField(atomic,ecm));
    }
}

void Nucleus::markWounded(Particle* part)
{
    if(part->getNumberOfCollision() == 1)
    {
        woundedNucleons.push_back(part);
    }
}

void Nucleus::clearNucleons()
{
    for(int i = 0; i < nucleons.size(); i++)
        delete nucleons[i];

    woundedNucleons.clear();
    nucleons.clear();
}

void Nucleus::GetDeuteronPosition(double& x1,double& y1,double& z1,double& x2,double& y2,double& z2)
{
   //get proton/neutron separation d (fm)
   double d;

   d=sample_deuteron.rand(); //now sample random rotation of deuteron
   x1 = (d/2.0);
   y1 = 0;
   z1 = 0;

   Point3D p3d(x1,y1,z1);
   p3d.rotate(ctr, phir);
   x1 = p3d.x; y1 = p3d.y; z1 = p3d.z;

   x2 = -x1;
   y2 = -y1;
   z2 = -z1;
}

void Nucleus::readin_triton_position()
{
   ifstream triton_position("tables/triton_positions.dat");
   double x1, y1, z1, x2, y2, z2, x3, y3, z3;
   triton_position >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
   while(!triton_position.eof())
   {
       vector<double> temp;
       temp.push_back(x1);
       temp.push_back(y1);
       temp.push_back(z1);
       temp.push_back(x2);
       temp.push_back(y2);
       temp.push_back(z2);
       temp.push_back(x3);
       temp.push_back(y3);
       temp.push_back(z3);
       triton_pos.push_back(temp);
       triton_position >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
   }
}

void Nucleus::readin_nucleon_positions()
{
   cout << "read in nucleon positions with nucleon-nucleon correlation..." << flush;
   ostringstream filename;
   if(atomic == 197)
   {
      filename << "tables/au197-sw-full_3Bchains-conf1820.dat";
      n_configuration = 1820;
   }
   else if(atomic == 208)
   {
      //int temp = rand() % 10 + 1;
      int temp = 1;
      filename << "tables/pb208-" << temp << ".dat";
      n_configuration = 10000;
   }
   else
      return;

   nucleon_pos_array = new double** [n_configuration];
   for(int iconf = 0; iconf < n_configuration; iconf++)
   {
      nucleon_pos_array[iconf] = new double* [atomic];
      for(int ia = 0; ia < atomic; ia++)
         nucleon_pos_array[iconf][ia] = new double [3];
   }

   ifstream input(filename.str().c_str());
   for(int iconf = 0; iconf < n_configuration; iconf++)
   {
      for(int ia = 0; ia < atomic; ia++)
      {
         double x_local, y_local, z_local;
         int isospin;
         int dummy;
         if(atomic == 208)
            input >> x_local >> y_local >> z_local >> isospin;
         else
            input >> x_local >> y_local >> z_local >> isospin >> dummy;
         nucleon_pos_array[iconf][ia][0] = x_local;
         nucleon_pos_array[iconf][ia][1] = y_local;
         nucleon_pos_array[iconf][ia][2] = z_local;
      }
   }
   input.close();
   cout << " done." << endl;
}

void Nucleus::GetTritonPosition(double& x1,double& y1,double& z1,double &x2,double& y2,double& z2, double &x3, double &y3, double &z3)
{
   int num_configuration = triton_pos.size();
   int rand_num = rand() % num_configuration;

   x1 = triton_pos[rand_num][0];
   y1 = triton_pos[rand_num][1];
   z1 = triton_pos[rand_num][2];
   x2 = triton_pos[rand_num][3];
   y2 = triton_pos[rand_num][4];
   z2 = triton_pos[rand_num][5];
   x3 = triton_pos[rand_num][6];
   y3 = triton_pos[rand_num][7];
   z3 = triton_pos[rand_num][8];

   Point3D p3d1(x1,y1,z1);
   p3d1.rotate(ctr, phir);
   x1 = p3d1.x; y1 = p3d1.y; z1 = p3d1.z;

   Point3D p3d2(x2,y2,z2);
   p3d2.rotate(ctr, phir);
   x2 = p3d2.x; y2 = p3d2.y; z2 = p3d2.z;

   Point3D p3d3(x3,y3,z3);
   p3d3.rotate(ctr, phir);
   x3 = p3d3.x; y3 = p3d3.y; z3 = p3d3.z;
}

// *** this function applies to deformed nuclei ***
void Nucleus::getDeformRandomWS(double& x, double& y, double& z)
{
  double rad1 = rad;
  double rwMax1 = rwMax;
  double r = 0.0;
  double cx = 1.0; // cx is cos(theta)

  if (deformed)
  {
    do {
      //Uniform distribution in a sphere with r = rmaxCut
      r = rmaxCut*pow(drand48(),1.0/3.0);
      cx = 1.0-2.0*drand48();

      //Deformation
      //Main axis in z-axis
      double y20 = SphericalHarmonics(2,cx);
      double y40 = SphericalHarmonics(4,cx);
      rad1 =rad*(1.0 + beta2 * y20 + beta4 * y40);
      rwMax1 =1.0/(1.0+exp(-rad1/dr));

    } while(drand48()*rwMax1 > 1.0/(1.0+exp((r-rad1)/dr)));

    double sx = sqrt(1.0-cx*cx); // sx is sin(theta)
    double phi = 2*M_PI*drand48();

    Point3D p3d(r*sx*cos(phi), r*sx*sin(phi), r*cx);
    p3d.rotate(ctr, phir);
    x = p3d.x; y = p3d.y; z = p3d.z;
  }
  else
  {
    do {
      //Uniform distribution in a sphere with r = rmaxCut
      r = rmaxCut*pow(drand48(),1.0/3.0);
    } while(drand48()*rwMax > 1.0/(1.0+exp((r-rad)/dr)));
    cx = 1.0-2.0*drand48();
    double sx = sqrt(1.0-cx*cx); // sx is sin(theta)
    double phi = 2*M_PI*drand48();
    x = r*sx*cos(phi);
    y = r*sx*sin(phi);
    z = r*cx;
  }
}

void Nucleus::get_nucleon_position_with_NN_correlation(double **nucleon_ptr)
{
   int i_configuration = rand() % n_configuration;  // pick the configuration
   double** temp_nucleus = new double* [atomic];
   for(int ia = 0; ia < atomic; ia++)
      temp_nucleus[ia] = new double [3];

   // first recenter the nucleus
   double xcm = 0.0, ycm = 0.0, zcm = 0.0;
   for(int ia = 0; ia < atomic; ia++)
   {
      xcm += nucleon_pos_array[i_configuration][ia][0];
      ycm += nucleon_pos_array[i_configuration][ia][1];
      zcm += nucleon_pos_array[i_configuration][ia][2];
   }
   for(int ia = 0; ia < atomic; ia++)
   {
      temp_nucleus[ia][0] = nucleon_pos_array[i_configuration][ia][0] - xcm/atomic;
      temp_nucleus[ia][1] = nucleon_pos_array[i_configuration][ia][1] - ycm/atomic;
      temp_nucleus[ia][2] = nucleon_pos_array[i_configuration][ia][2] - zcm/atomic;
   }

   // then rotate a random angle
   double cos_theta = 1.0-2.0*drand48();
   double phi = 2*M_PI*drand48();
   for(int ia = 0; ia < atomic; ia++)
   {
      Point3D p3d(temp_nucleus[ia][0], temp_nucleus[ia][1], temp_nucleus[ia][2]);
      p3d.rotate(cos_theta, phi);
      temp_nucleus[ia][0] = p3d.x; 
      temp_nucleus[ia][1] = p3d.y; 
      temp_nucleus[ia][2] = p3d.z;
   }

   // return back
   for(int ia = 0; ia < atomic; ia++)
      for(int ii = 0; ii < 3; ii++)
         nucleon_ptr[ia][ii] = temp_nucleus[ia][ii];

   // clean up
   for(int ia = 0; ia < atomic; ia++)
      delete [] temp_nucleus[ia];
   delete [] temp_nucleus;
}


double Nucleus::SphericalHarmonics(int l, double ct)
{
  //Currently assuming m=0 and available for Y_{20} and Y_{40}
  // "ct" is cos(theta)

  double ylm = 0.0;

  if(l == 2){

    ylm = 3.0*ct*ct-1.0;
    ylm *= 0.31539156525252005; //pow(5.0/16.0/M_PI,0.5);

  }else if (l == 4){

    ylm = 35.0*ct*ct*ct*ct;
    ylm -= 30.0*ct*ct;
    ylm += 3.0;
    ylm *= 0.10578554691520431; //3.0/16.0/pow(M_PI,0.5);

  }else{
    cerr << "Not available in Overlap::SphericalHarmonics" << endl;
  }

  return ylm;

}

//**********************************************************************
void Nucleus::Gauss38(double xini,double xfin,double* xn,double* wn)
{
    double  x[38], w[38];

    x[37]=9.980499305357e-1;
    x[36]=9.897394542664e-1;
    x[35]=9.748463285902e-1;
    x[34]=9.534663309335e-1;
    x[33]=9.257413320486e-1;
    x[32]=8.918557390046e-1;
    x[31]=8.520350219324e-1;
    x[30]=8.065441676053e-1;
    x[29]=7.556859037540e-1;
    x[28]=6.997986803792e-1;
    x[27]=6.392544158297e-1;
    x[26]=5.744560210478e-1;
    x[25]=5.058347179279e-1;
    x[24]=4.338471694324e-1;
    x[23]=3.589724404794e-1;
    x[22]=2.817088097902e-1;
    x[21]=2.025704538921e-1;
    x[20]=1.220840253379e-1;
    x[19]=4.078514790458e-2;

//    .....   WEIGHT       ...........
    w[37]=5.002880749632e-3;
    w[36]=1.161344471647e-2;
    w[35]=1.815657770961e-2;
    w[34]=2.457973973823e-2;
    w[33]=3.083950054518e-2;
    w[32]=3.689408159400e-2;
    w[31]=4.270315850467e-2;
    w[30]=4.822806186076e-2;
    w[29]=5.343201991033e-2;
    w[28]=5.828039914700e-2;
    w[27]=6.274093339213e-2;
    w[26]=6.678393797914e-2;
    w[25]=7.038250706690e-2;
    w[24]=7.351269258474e-2;
    w[23]=7.615366354845e-2;
    w[22]=7.828784465821e-2;
    w[21]=7.990103324353e-2;
    w[20]=8.098249377060e-2;
    w[19]=8.152502928039e-2;

    for(int i=0;i<19;i++) {
      x[i] = -x[37-i];
      w[i] =  w[37-i];
    }
    for(int i=0;i<38;i++) {
      xn[i] =(xfin-xini)*x[i]/2.0+(xini+xfin)/2.0;
      wn[i] =(xfin-xini)*w[i]/2.0;
    }

}

void Nucleus::dumpParticipants(ofstream &of)
{
    double x, y;
    for(int i = 0; i < woundedNucleons.size(); i++)
    {
        x = woundedNucleons[i]->getX();
        y = woundedNucleons[i]->getY();
        of  << setprecision(3) << setw(10) << x << "   "
            << setprecision(3) << setw(10) << y << "   " 
            << nuclID << endl;
    }
}

void Nucleus::dumpNucleons(ofstream& of)
{
    
    double x, y;
    for(int i = 0; i < nucleons.size(); i++)
    {
        x = nucleons[i]->getX();
        y = nucleons[i]->getY();
        of  << setprecision(3) << setw(10) << x
            << setprecision(3) << setw(10) << y
            << endl;
    }
}

void Nucleus::dumpQuarks(ofstream& of)
{
    double x, y;
    vector<Quark> quarks;
    for(int i = 0; i < woundedNucleons.size(); i++)
    {
        quarks = woundedNucleons[i]->getQuarks();
        for(int j = 0; j < quarks.size(); j++)
        {
            x = quarks[j].getX();
            y = quarks[j].getY();
            of << setprecision(3) << setw(10) << x
               << setprecision(3) << setw(10) << y << " "
               << quarks[j].getBoundingBox().toString() << endl;
               
        }
    }
}



