#include "MCnucl.h"
#include "Regge96.h"
#include <iomanip>
#include <ctime>
#include <sys/time.h>
#include <memory>

#include "arsenal.h"
#include "GaussianNucleonsCal.h"
#include "ParameterReader.h"
#include "NBD.h"

#define HBARC 0.197327053

using namespace std;

//*********************************************************************
// functions to initialize the 2d transv. coordinate grid,
// place two colliding nuclei on the grid, determine their thickness
// generate lookup table for dN/dy as fct of T_A(rt), T_B(rt)

MCnucl::MCnucl(ParameterReader* paraRdr_in)
{
  paraRdr = paraRdr_in;

  // default Alpha is Glauber model: sd = (1-Alpha)WN + (Alpha)BC
  Alpha = paraRdr->getVal("alpha");

  // tmax-1 is max # of locally overlapping nucleons, upper limit in dNdy table
  tmax = paraRdr->getVal("tmax");

  // fix grid properties
  Xmax = paraRdr->getVal("maxx");
  Ymax = paraRdr->getVal("maxy");
  Xmin = -Xmax;
  Ymin = -Ymax;
  dx = paraRdr->getVal("dx");
  dy = paraRdr->getVal("dy");
  Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
  Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;

  // default npart cutoff
  NpartMax = paraRdr->getVal("Npmax");
  NpartMin = paraRdr->getVal("Npmin");

  // Unintegrated PT optns
  PTinte = paraRdr->getVal("PT_flag");
  PTmax  = paraRdr->getVal("PT_Max");
  PTmin  = paraRdr->getVal("PT_Min");
  dpt = paraRdr->getVal("d_PT");
  MaxPT=(int)((PTmax-PTmin)/dpt+0.1)+1;
  if(PTinte<0)
      PT_order = paraRdr->getVal("PT_order");   
  else
      PT_order = 1; //does not apply when there is no PT integration
  
  //.... NN cross sections in mb
  double ecm = paraRdr->getVal("ecm");
  double sig = hadronxsec::totalXsection(200.0,0);
  double sigel = hadronxsec::elasticXsection(sig,200.0,0,0);
  siginNN200 = sig - sigel;
  sig = hadronxsec::totalXsection(ecm,0);
  sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
  siginNN = sig - sigel;

  // Include additional CC fluctution?
  CCFluctuationModel = paraRdr->getVal("cc_fluctuation_model");
  CCFluctuationK = paraRdr->getVal("cc_fluctuation_k");
  if (CCFluctuationModel) nbd = new NBD;
  if (CCFluctuationModel > 5)
  {
     gsl_rng_env_setup();
     gslRngType = gsl_rng_default;
     gslRng = gsl_rng_alloc(gslRngType);
     timeval a;
     gettimeofday(&a, 0);
     int randomSeed=a.tv_usec; // randomSeed use CPU clock
     gsl_rng_set (gslRng, (unsigned long int) randomSeed); //initialize random generator
     ccFluctuationGammaTheta = paraRdr->getVal("cc_fluctuation_Gamma_theta");
  }


  which_mc_model = paraRdr->getVal("which_mc_model");
  sub_model = paraRdr->getVal("sub_model");
  shape_of_nucleons = paraRdr->getVal("shape_of_nucleons");
  

  gaussCal = NULL;
  entropy_gaussian_width = 0.0;
  paraRdr->setVal("siginNN", siginNN);
  gaussCal = new GaussianNucleonsCal(paraRdr); // for Gaussian-shaped nucleons calculations
  entropy_gaussian_width = gaussCal->width;
  entropy_gaussian_width_sq = entropy_gaussian_width*entropy_gaussian_width;
  
  
  proj = new Nucleus(paraRdr->getVal("Aproj"),
                    paraRdr,
                    paraRdr->getVal("proj_deformed"));
  targ = new Nucleus(paraRdr->getVal("Atarg"),
                    paraRdr,
                    paraRdr->getVal("targ_deformed"));

  // adding quark substructure Fluctuations (from Kevin Welsh)
  shape_of_entropy = paraRdr->getVal("shape_of_entropy");
  forceCollisionCriterion = paraRdr->getVal("collision_criterion");

  dndyTable=0;    // lookup table pointers not valid yet
  dndydptTable=0;
  overSample=1;  // default: no oversampling
  binRapidity = paraRdr->getVal("ny");
  rapMin = paraRdr->getVal("rapMin");
  rapMax = paraRdr->getVal("rapMax");
  rapidity = 0.0;

  dsq = 0.1*siginNN/M_PI/overSample;

  kln=0;  // pointer to class for small-x gluons
  val=0;  // pointer to class for large-x (x>x0) partons

  TA1 = new double* [Maxx];    // 2d grid for proj/targ. thickness functions
  TA2 = new double* [Maxx];
  rho_binary = new double* [Maxx];  // 2d grid for the binary density
  for(int ix=0;ix<Maxx;ix++)
  {
    TA1[ix] = new double [Maxy];
    TA2[ix] = new double [Maxy];
    rho_binary[ix] = new double [Maxy];
    for(int iy = 0; iy < Maxy; iy++)
    {
       TA1[ix][iy] = 0;
       TA2[ix][iy] = 0;
       rho_binary[ix][iy] = 0;
    }
  }

  rho = new GlueDensity(Xmax,Ymax,PTmin,PTmax,dx,dy,dpt,binRapidity,rapMin,rapMax);

  Xcm=0.0, Ycm=0.0, angPart=0.0;

  // flag to include nucleon-nucleon correlation when sampled the nucleon positions
  flag_NN_correlation = paraRdr->getVal("include_NN_correlation");

}


MCnucl::~MCnucl()

{
  delete proj;
  delete targ;
  
  for(int ix = 0; ix < Maxx; ix++)
  {
    delete [] TA1[ix];
    delete [] TA2[ix];
    delete [] rho_binary[ix];
  }
  delete [] TA1;
  delete [] TA2;
  delete [] rho_binary;

  delete  rho;

  if(dndyTable) {
    for(int iy=0;iy<binRapidity;iy++) {
      for(int j=0;j<tmax;j++) delete [] dndyTable[iy][j];
      delete [] dndyTable[iy];
    }
    delete [] dndyTable;
  }

  if(dndydptTable) {
    for (int iy=0; iy<binRapidity; iy++) {
      for (int j=0; j<tmaxPt; j++) {
        for (int i=0; i<tmaxPt; i++) delete [] dndydptTable[iy][j][i];
        delete [] dndydptTable[iy][j];
      }
      delete [] dndydptTable[iy];
    }
    delete [] dndydptTable;
  }

  if(CCFluctuationModel > 5) gsl_rng_free(gslRng);
  if(nbd) delete nbd;
  if (gaussCal) delete gaussCal;
}



/* place two nuclei/nucleons on the transverse lattice (separated by
   impact parameter b)
*/
void MCnucl::generateNuclei(double b)
{
    for(int ie=0;ie<overSample;ie++) {
        proj->populate(b/2.0,0);
        targ->populate(-b/2.0,0);
    }
}

// --- find participants from proj/target and the number of binary coll. ---
int MCnucl::getBinaryCollision()
{
  bool missingNucleus = false;
  
  // Handling for the intrinsic nucleus case
  if(proj->getAtomic() == 0)
  {
      vector<Particle*> nucl2 = targ->getNucleons();
      missingNucleus = true;
      cout << "Projectile missing. Dumping target entropy." << endl;
      for(int i = 0; i < (int)nucl2.size(); i++){
          selectFluctFactors(nucl2[i]);
          targ->markWounded(nucl2[i]);
      }
  }
  else if(targ->getAtomic() == 0)
  {
      vector<Particle*> nucl1 = proj->getNucleons();
      missingNucleus = true;
      cout << "Target missing. Dumping projectile entropy." << endl;
      for(int i = 0; i < (int)nucl1.size(); i++){
          selectFluctFactors(nucl1[i]);
          proj->markWounded(nucl1[i]);
      }
  }
  else
  {
      vector<BoundingBox> projBoxes = proj->getSortedBoundingBoxes();
      vector<BoundingBox> targBoxes = targ->getSortedBoundingBoxes();
      
      int startingIndex = 0;
      for(int iproj = 0; iproj < projBoxes.size(); iproj++)
      {
          BoundingBox projBox = projBoxes[iproj];
          BoundingBox targBox = targBoxes[startingIndex];
          
          // Skip the left most nucleons for each proj nucleon.
          while(targBox.getXR() < projBox.getXL()
                  && startingIndex < targBoxes.size())
          {
              startingIndex++;
              targBox = targBoxes[startingIndex];
          }
          
          // Actually test a collision for the next ones,
          // until they get too far away.
          int i = startingIndex;
          while(projBox.getXR() >= targBox.getXL()
                  && i < targBoxes.size())
          {
              if(projBox.getYL() <= targBox.getYR())
              {
                  if(projBox.getYR() >= targBox.getYL())
                  {
                      // Now we know the boxes do overlap in x and y.
                      Particle* projPart = projBox.getParticle();
                      Particle* targPart = targBox.getParticle();
                      if(hit(projPart,targPart,projBox.intersection(targBox)))
                      {
                          selectFluctFactors(projPart);
                          selectFluctFactors(targPart);
                          projPart->addCollidingParticle(targPart);
                          targPart->addCollidingParticle(projPart);
                          proj->markWounded(projPart);
                          targ->markWounded(targPart);
                      }
                  }
              }
              i++;
              targBox = targBoxes[i];
          }
          
          // Continue loop over projBoxes
      }
      
      // Exit hit detection
  }

  createBinaryCollisions();
  Npart1=proj->getNpart();
  Npart2=targ->getNpart();
  
  Npart1 /= overSample;
  Npart2 /= overSample;
  
  return missingNucleus ? 1 : binaryCollision.size();
}

void MCnucl::selectFluctFactors(Particle* part)
{
    if(CCFluctuationModel > 5)
    {
      if(shape_of_entropy == 3)
      {
        double f1 = sampleFluctuationFactorforParticipant();
        double f2 = sampleFluctuationFactorforParticipant();
        double f3 = sampleFluctuationFactorforParticipant();
        part->setQuarkFluctfactor(f1,f2,f3);
      }
      else
        part->setFluctfactor(sampleFluctuationFactorforParticipant());
    }
}

void MCnucl::createBinaryCollisions()
{
    vector<Particle*> projParticipants = proj->getParticipants();
    // Loop through the participants for just the projectile,
    // then loop through all the particles that collided with each.
    for(int i = 0; i < projParticipants.size(); i++)
    {
        Particle* part = projParticipants[i];
        double partX = part->x;
        double partY = part->y;
        
        vector<Particle*> collidingParticles = part->getCollidingParticles();
        for(int j = 0; j < collidingParticles.size(); j++)
        {
            Particle* colliding = collidingParticles[j];
            CollisionPair* pair = 
                    new CollisionPair((partX+colliding->x)/2.0,
                                      (partY+colliding->y)/2.0);
            
            if(CCFluctuationModel > 5)
                pair->setfluctfactor(sampleFluctuationFactorforBinaryCollision());
            
            if(which_mc_model == 5 && sub_model == 2){
                pair->additional_weight += 1/part->getNumberOfCollision();
                pair->additional_weight += 1/colliding->getNumberOfCollision();
            }
            binaryCollision.push_back(pair);
        }
    }
}

// old stuff
//  (YN): determine whether nucleons separated by distance dr2 interact or not
// interaction probability at impact param. b is 1-exp(-sigeff(s)*Tpp(b))
int MCnucl::hit(Particle* part1, Particle* part2,const Box2D& overlapRegion)
{
    double b = sqrt((part2->x-part1->x)*(part2->x-part1->x)
                +(part2->y-part1->y)*(part2->y-part1->y));
    
    switch(forceCollisionCriterion)
    {
        case 1:
            disk:
            return (b*b<=dsq) ? 1 : 0;
        case 2:
            gauss:
            return gaussCal->testSmoothCollision(b);
        default:
            if(shape_of_nucleons == 1)
                goto disk;
            if(shape_of_entropy == 2)
                goto gauss;
            return gaussCal->testFluctuatedCollision(part1,part2,overlapRegion);
    }
}
// checks whether Npart1+Npart2 is in the desired range
//   ATTN: Npart1, Npart2 need to be initialized through getBinaryCollision() !
int MCnucl::CentralityCut()
{
  int Nptot = Npart1 + Npart2;
  if (Nptot<=NpartMax && Nptot>NpartMin) return 1;
  return 0;
}


// --- determine thickness of proj+targ nuclei over 2d transv. grid,
//     for given MC event ---
void MCnucl::getTA2()
{
  int imax=0;
  
  for(int ix=0;ix<Maxx;ix++){
      for(int iy=0;iy<Maxy;iy++)
      {
          TA1[ix][iy]=0.0;
          TA2[ix][iy]=0.0;
      }
  }
  
  setThickness(proj,TA1);
  setThickness(targ,TA2);

  for(int ix=0;ix<Maxx;ix++)
      for(int iy=0;iy<Maxy;iy++)
      {
          int i = int((TA1[ix][iy])/(10.0/siginNN)+0.5);  // keep track of highest density
          int j = int((TA2[ix][iy])/(10.0/siginNN)+0.5);
          imax = max(imax,i);
          imax = max(imax,j);
      }

  if(which_mc_model == 1)
  {
     if(imax>=tmax) {
       cout  << "# WARNING: in MCnucl::getTA2() : imax=" << imax
             << " should be less than tmax=" << tmax << endl;
       exit(0);
     }
  }
}

void MCnucl::setThickness(Nucleus* nucl, double ** TA)
{
    vector<Particle*> participant = nucl->getParticipants();
    
    double nucleon_width;
    if (shape_of_nucleons>=2 && shape_of_nucleons<=9) 
        nucleon_width = gaussCal->width;
    
    double d_max;
    if (shape_of_nucleons == 1)
        d_max = 2.*sqrt(dsq);
    if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
        d_max = 5.*nucleon_width;
    
    for(unsigned int ipart=0; ipart<participant.size(); ipart++)
    {
      double x = participant[ipart]->x;
      double y = participant[ipart]->y;
      int x_idx_left = (int)((x - d_max - Xmin)/dx);
      int x_idx_right = (int)((x + d_max - Xmin)/dx);
      int y_idx_left = (int)((y - d_max - Ymin)/dy);
      int y_idx_right = (int)((y + d_max - Ymin)/dy);
      x_idx_left = max(0, x_idx_left);
      x_idx_right = min(Maxx, x_idx_right);
      y_idx_left = max(0, y_idx_left);
      y_idx_right = min(Maxy, y_idx_right);
      for(int ix = x_idx_left; ix < x_idx_right; ix++)
      {
         double xg = Xmin + ix*dx;
         for(int iy = y_idx_left; iy < y_idx_right; iy++)
         {
             double yg = Ymin + iy*dy;
             double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
             if (shape_of_nucleons==1) // "Checker" nucleons:
             {
               if(dc>dsq) 
                   continue;
               TA[ix][iy] += 10.0/siginNN;
             }
             else if (shape_of_nucleons>=2 && shape_of_nucleons<=9) // Gaussian nucleons:
             {
                TA[ix][iy] += participant[ipart]->getSmoothTn(xg,yg);
             }
         }
      }
    }
}

// calculate the binary collision density in the transverse plane
void MCnucl::calculate_rho_binary()
{
    // We only need to look at one of the nuclei's wounded nucleons
    // since all colliding nucleons in the targ are in who_hit_me for the
    // proj.
   double dc_sq_max_gaussian = 25.*entropy_gaussian_width_sq;
   double d_max;
   if (shape_of_nucleons == 1)
       d_max = 2.*sqrt(dsq);
   if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
       d_max = 5.*entropy_gaussian_width;
   
   for(int ir = 0; ir < Maxx; ir++)
       for(int jr = 0; jr < Maxy; jr++)
           rho_binary[ir][jr] = 0.0;
   
   for(int i = 0; i < binaryCollision.size(); i++)
   {
       double x = binaryCollision[i]->x;
       double y = binaryCollision[i]->y;
       int x_idx_left = (int)((x - d_max - Xmin)/dx);
       int x_idx_right = (int)((x + d_max - Xmin)/dx);
       int y_idx_left = (int)((y - d_max - Ymin)/dy);
       int y_idx_right = (int)((y + d_max - Ymin)/dy);
       x_idx_left = max(0, x_idx_left);
       x_idx_right = min(Maxx, x_idx_right);
       y_idx_left = max(0, y_idx_left);
       y_idx_right = min(Maxy, y_idx_right);
       for(int ir = x_idx_left; ir < x_idx_right; ir++)
       {
          double xg = Xmin + ir*dx;
          for(int jr = y_idx_left; jr < y_idx_right; jr++)
          {
             double yg = Ymin + jr*dy;
             double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
             if (shape_of_nucleons == 1)
             {
                if(dc <= dsq) 
                   rho_binary[ir][jr] += (10.0/siginNN);
             }
             else if (shape_of_nucleons>=2 && shape_of_nucleons<=9)
             {
                if (dc > dc_sq_max_gaussian) 
                   continue; // skip small numbers to speed up
                rho_binary[ir][jr] += GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2*entropy_gaussian_width_sq)); 
                // this density is normalized to 1, to be consistent with the disk-like treatment; 
             }
          }
       }
   }
}

// --- initializes dN/dyd2rt (or dEt/...) on 2d grid for rapidity slice iy
//     and integrates it to obtain dN/dy (or dEt/dy) ---
void MCnucl::setDensity(int iy, int ipt)
{
  // which_mc_model==1 -> KLN-like
  if (which_mc_model==1 && ipt>=0 && (dndydptTable==0)) {
    cout << "ERROR in MCnucl::setDensity() : pt-bin=" << ipt <<
      " but no dndydptTable !" << endl;
    exit(0);
  }

  // which_mc_model==1 -> KLN-like
  if (which_mc_model==1 && ipt<0 && (dndyTable==0)) {
    cout <<
     "ERROR in MCnucl::setDensity() : pt-integrated yields require dndyTable !" << endl;
    exit(0);
  }

  double tblmax=0, table_result=0;

  rapidity=rapMin + (rapMax-rapMin)/binRapidity*iy;
  dndy=0.0;
 
  double dc_sq_max_gaussian = 25.*entropy_gaussian_width_sq;
  double d_max;
  if (shape_of_nucleons == 1)
      d_max = 2.*sqrt(dsq);
  if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
      d_max = 5.*entropy_gaussian_width;
  
  if(which_mc_model==1) // MC-KLN
  {
    for(int ir=0;ir<Maxx;ir++)  // loop over 2d transv. grid
      for(int jr=0;jr<Maxy;jr++) 
      {
        double xg = Xmin + ir*dx;
        double yg = Ymin + jr*dy;
        double di = TA1[ir][jr]/dT;
        double dj = TA2[ir][jr]/dT;
        tblmax=max(di,tblmax);
        tblmax=max(dj,tblmax);
        if((di<0 || di>=tmax-2) || (dj<0 || dj>=tmax-2) ) {
            cerr << "di= " << di << " dj= " << dj
                << " You should increase the dimension of dndyTable"
                << endl;
            exit(1);
        }
        int i = floor(di); int j = floor(dj);
        if (ipt<0) // without pt dependence
        {
          table_result = sixPoint2dInterp(di-i, dj-j, // x and y value, in lattice unit (dndyTable_step -> 1)
          dndyTable[iy][i][j], dndyTable[iy][i][j+1], dndyTable[iy][i][j+2], dndyTable[iy][i+1][j], dndyTable[iy][i+1][j+1], dndyTable[iy][i+2][j]);
          rho->setDensity(iy,ir,jr,table_result);
        }
        else // with pt dependence
        {
          table_result = sixPoint2dInterp(di-i, dj-j, // x and y value, in lattice unit (dndyTable_step -> 1)
          dndydptTable[iy][i][j][ipt], dndydptTable[iy][i][j+1][ipt], dndydptTable[iy][i][j+2][ipt], dndydptTable[iy][i+1][j][ipt], dndydptTable[iy][i+1][j+1][ipt], dndydptTable[iy][i+2][j][ipt]);
          rho->setDensity(iy,ir,jr, ipt, table_result);
        }
        //dndy += dndyTable[iy][i][j];
        dndy += table_result;
      }
  }
  else if(which_mc_model==5) // MC-Glb, entropy centered at wounded nucleons and binary collision centers
  {
      double **rhop, **tab;
      rhop = new double * [Maxx];
      tab = new double * [Maxx];
      for(int ir = 0; ir < Maxx; ir++)
      {
          rhop[ir] = new double [Maxy];
          tab[ir] = new double [Maxy];
          for(int jr = 0; jr < Maxy; jr++)
          {
              rhop[ir][jr] = 0.0;
              tab[ir][jr] = 0.0;
          }
      }
      // wounded nucleon treatment:
      if (sub_model==1) // "classical" Glb
      {
          addEntropyDensity(proj,rhop);
          addEntropyDensity(proj,rhop);
          //rhop = (TA1[ir][jr]+TA2[ir][jr])*(1.0-Alpha)/2;
          
          double prefactor = (1.0 - Alpha)/2.;
          for(int ir = 0; ir < Maxx; ir++)
              for(int jr = 0; jr < Maxy; jr++)
                  rhop[ir][jr] = rhop[ir][jr]*prefactor;
      }
      else if (sub_model==2) // "Ulrich" Glb
      {
          for(int ir = 0; ir < Maxx; ir++)
              for(int jr = 0; jr < Maxy; jr++)
                  rhop[ir][jr] = 0.0; // no "wounded" contribution
      }
      else
      {
          cout << "MCnucl::setDensity error: which_mc_model is set to " << which_mc_model << ", but the associated sub_model " << sub_model << " is not recognized." << endl;
          exit(-1);
      }

      // binary collision treatment:
      if(Alpha > 1e-8)
      {
          double fluctfactor = 1.0;
          for(int icoll=0;icoll<binaryCollision.size();icoll++)
          {
              double x = binaryCollision[icoll]->x;
              double y = binaryCollision[icoll]->y;
              if(CCFluctuationModel > 5)
                  fluctfactor = binaryCollision[icoll]->getfluctfactor();
              int x_idx_left = (int)((x - d_max - Xmin)/dx);
              int x_idx_right = (int)((x + d_max - Xmin)/dx);
              int y_idx_left = (int)((y - d_max - Ymin)/dy);
              int y_idx_right = (int)((y + d_max - Ymin)/dy);
              x_idx_left = max(0, x_idx_left);
              x_idx_right = min(Maxx, x_idx_right);
              y_idx_left = max(0, y_idx_left);
              y_idx_right = min(Maxy, y_idx_right);
              for(int ir = x_idx_left; ir < x_idx_right; ir++)
              {
                 double xg = Xmin + ir*dx;
                 for(int jr = y_idx_left; jr < y_idx_right; jr++)
                 {
                     double yg = Ymin + jr*dy;
                     double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
                     if (shape_of_entropy==1)
                     {
                         if(dc <= dsq) 
                             tab[ir][jr] += fluctfactor*(10.0/siginNN)*(Alpha + (1.-Alpha)*binaryCollision[icoll]->additional_weight); // second part in the paranthesis is for Uli-Glb model
                     }
                     else if (shape_of_entropy>=2 && shape_of_entropy<=9)
                     {
                         if (dc <= dc_sq_max_gaussian)
                             tab[ir][jr] += fluctfactor*GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2*entropy_gaussian_width_sq))*(Alpha + (1.-Alpha)*binaryCollision[icoll]->additional_weight); // this density is normalized to 1, to be consisitent with the disk-like treatment; second part in the last parathesis is for Uli-Glb model
                     }
                 }
              }
          }
      } 
      else
      {
          for(int ir = 0; ir < Maxx; ir++)
              for(int jr = 0; jr < Maxy; jr++)
                  tab[ir][jr] = 0.;
      } // if(Alpha > 1e-8)

      // set entropy density
      for(int ir = 0; ir < Maxx; ir++)
          for(int jr = 0; jr < Maxy; jr++)
          {
              double density = rhop[ir][jr] + tab[ir][jr]; // "rhop" (n_WN) and "tab" (n_bin) have already been multiplied by (1-x)/2 and x (x=Alpha here) correspondingly
              rho->setDensity(iy,ir,jr,density);
              dndy += density;
          }
      for(int ir = 0; ir < Maxx; ir++)
      {
          delete [] rhop[ir];
          delete [] tab[ir];
      }
      delete [] rhop;
      delete [] tab;
  }
  if(dndy<1e-15)
  {
    cout << "MCnucl::setDensity dndy = 0 !!  y= " << rapidity
     << " dndy= " << dndy << endl;
    exit(0);
  }

  // Should I include additional fluctuation for MCKLN?
  if (CCFluctuationModel>0 && CCFluctuationModel <= 5) fluctuateCurrentDensity(iy);
}

void MCnucl::addEntropyDensity(Nucleus* nucl, double** dens)
{
    vector<Particle*> participant = nucl->getParticipants();
    double d_max;
    if (shape_of_nucleons == 1)
        d_max = 2.*sqrt(dsq);
    if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
        d_max = 5.*entropy_gaussian_width;
    for(unsigned int ipart=0; ipart<participant.size(); ipart++) 
    {
      double x = participant[ipart]->x;
      double y = participant[ipart]->y;
      int x_idx_left = (int)((x - d_max - Xmin)/dx);
      int x_idx_right = (int)((x + d_max - Xmin)/dx);
      int y_idx_left = (int)((y - d_max - Ymin)/dy);
      int y_idx_right = (int)((y + d_max - Ymin)/dy);
      x_idx_left = max(0, x_idx_left);
      x_idx_right = min(Maxx, x_idx_right);
      y_idx_left = max(0, y_idx_left);
      y_idx_right = min(Maxy, y_idx_right);


      for(int ir = x_idx_left; ir < x_idx_right; ir++)
      {
         double xg = Xmin + ir*dx;
         for(int jr = y_idx_left; jr < y_idx_right; jr++)
         {
             double yg = Ymin + jr*dy;
             double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
             if (shape_of_entropy==1) // "Checker" nucleons:
             {
               if(dc>dsq) 
                   continue;

               double areai = 10.0/siginNN;
               dens[ir][jr] += areai*participant[ipart]->getFluctfactor();
             }
             else if (shape_of_entropy>=2 && shape_of_entropy<=9) // Gaussian nucleons:
             {
               double density;
               if (shape_of_entropy == 3)
               {
                  density = participant[ipart]->getFluctuatedDensity(xg,yg);
               }
               else
               {
                  density = participant[ipart]->getSmoothDensity(xg,yg);
               }
               dens[ir][jr] += density;
             }
         }
      }
    }
}

//----------------------------------------------------------------------
void MCnucl::fluctuateCurrentDensity(int iy)
// Fluctuate the density profile in rho
{
    if (CCFluctuationModel==1) // use constant k/n:
    {
        for(int ir=0;ir<Maxx;ir++)  // loop over 2d transv. grid
        for(int jr=0;jr<Maxy;jr++)
        {
            double nb = rho->getDensity(iy,ir,jr)*dx*dy;
            double n = nbd->rand(nb/(nb+CCFluctuationK), CCFluctuationK);
            rho -> setDensity(iy,ir,jr,n/(dx*dy));
        }
    }
    else if (CCFluctuationModel==2)
    {
        // double kpp = 1.0/M_PI*dx*dy*1.0*kln->getLambdaQCD(); // second 1.0 is delta-eta
        double kpp = 1.0/M_PI*dx*dy*1.0*(0.25*0.25/HBARC/HBARC); // lambdaQCD=0.25/hbarC
        for(int ir=0;ir<Maxx;ir++)  // loop over 2d transv. grid
        for(int jr=0;jr<Maxy;jr++)
        {
            double k = kpp*min(TA1[ir][jr], TA2[ir][jr])*siginNN/10;
            double nb = rho->getDensity(iy,ir,jr)*dx*dy;
            double n;
            if (nb<1e-10)
                n = nb;
            else
                n = nbd->rand(nb/(nb+k), k);
            rho -> setDensity(iy,ir,jr,n/(dx*dy));
        }
    }
    else
    {
        cout << "MCnucl::fluctuateCurrentDensity error: CCFluctuationModel not supported." << endl;
        cout << "MCnucl:: CCFluctuationModel = " << CCFluctuationModel << endl;
        exit(-1);
    }
}


/***
   generates lookup table for dN/dy for varying proj/target thicknesses
***/
void MCnucl::makeTable()
{
  cout << "MCnucl::makeTable(): precalculating dNdy for all combinations of Ta and Tb." << endl;

  dT=10.0/siginNN/overSample;   // elementary thickness step
  if (shape_of_nucleons>=2 && shape_of_nucleons<=9) {  // Gaussian nucleons require finer steps
    tmax = paraRdr->getVal("tmax_subdivision")*(tmax-1) + 1;
    dT /= paraRdr->getVal("tmax_subdivision");
  }

  // allocate table
  dndyTable = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dndyTable[iy] = new double* [tmax];
    for(int j=0;j<tmax;j++) dndyTable[iy][j] = new double [tmax];
  }

int progress_counter = 0, progress_percent = 0, last_update = 0;
//===========================================================================
  for(int iy=0;iy<binRapidity;iy++) { // loop over rapidity bins
    double y = rapMin+(rapMax-rapMin)/binRapidity*iy;

    for(int i=0;i<tmax;i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for(int j=0;j<tmax;j++) { // loop over targ thickness
        double ta2 = dT*j;
        if(i>0 && j>0) {  // store corresponding dN/dy in lookup table
          // small-x gluons via kt-factorization
          dndyTable[iy][i][j] = kln->getdNdy(y,ta1,ta2, -1, PT_order); 
          // add large-x partons via DHJ formula if required
          if (val)
            dndyTable[iy][i][j] += val->getdNdy(y,ta1,ta2);
          //cout << ta1 << ", " << ta2 << ", " << dndyTable[iy][i][j] << endl;
        } else dndyTable[iy][i][j] = 0.0;
      progress_counter++;
      progress_percent = (progress_counter*100) / (binRapidity*tmax*tmax);
      if(((progress_percent%10) == 0) && (progress_percent != last_update))
        {
       cout << progress_percent << "% : " << std::flush;
       last_update = progress_percent;
        }
      }
    }
  }
  cout << endl;
  //===========================================================================
  cout << "MCnucl::makeTable(): done" << endl << endl;

  dumpdNdyTable4Col("data/dNdyTable.dat", dndyTable, 0);
}



/***
   generates lookup table for dN/dyd2pt for varying proj/target thicknesses
***/
void MCnucl::makeTable(double ptmin, double dpt, int iPtmax)
{
  cout << "MCnucl::makeTable(double, double, int): generating dN/dyd2pt lookup table... " << endl;

  dT=10.0/siginNN/overSample;   // elementary thickness step
  if (shape_of_nucleons>=2 && shape_of_nucleons<=9) {  // Gaussian nucleons require finer steps
    tmax = paraRdr->getVal("tmax_subdivition")*(tmax -1 ) + 1;
    dT /= paraRdr->getVal("tmax_subdivition");
  }

  // range of thicknesses for pt dependent lookup table
  tmaxPt = tmax;
  //tmaxPt = tmax/5;   // sufficient for pp

  // allocate table
  iptmax = iPtmax;  // destructor needs to delete table, store size
  dndydptTable = new double*** [binRapidity];
  for (int iy=0; iy<binRapidity; iy++) {
    dndydptTable[iy] = new double** [tmaxPt];
    for (int j=0; j<tmaxPt; j++) {
      dndydptTable[iy][j] = new double* [tmaxPt];
      for (int i=0; i<tmaxPt; i++)
        dndydptTable[iy][j][i] = new double [iptmax];
    }
  }

  
  int progress_counter = 0, progress_percent = 0, last_update = 0;
  //===========================================================================

  for(int iy=0;iy<binRapidity;iy++) { // loop over rapidity bins
    double y = rapMin+(rapMax-rapMin)/binRapidity*iy;
    for (int i=0; i<tmaxPt; i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for (int j=0; j<tmaxPt; j++) { // loop over targ thickness
        double ta2 = dT*j;
        if (val) val->lgXClearPtTable(); // recomp. large-rap pt-distributions
        for (int ipt=0; ipt<iptmax; ipt++) { // loop over pt bins
          if(i>0 && j>0) {  // store corresponding dN/dyd2pt in lookup table
            // high-rap *hadrons* via DHJ formula
            if (val)  dndydptTable[iy][i][j][ipt] =
            val->getdNdyd2pt(y,ta1,ta2,ptmin+ipt*dpt);
            else
              // small-x gluons via kt-factorization;  fixed pt, no integration
              dndydptTable[iy][i][j][ipt] = kln->getdNdy(y,ta1,ta2,ptmin+ipt*dpt);
          } else dndydptTable[iy][i][j][ipt] = 0.0;        
          progress_counter++;
          progress_percent = (progress_counter*100) / (binRapidity*tmax*tmax*iptmax);
          if(((progress_percent%10) == 0) && (progress_percent != last_update))
          {
           cout << progress_percent << "% : " << std::flush;
           last_update = progress_percent;
          }
        }
      }
    }
  }
  cout << "MCnucl::makeTable(double, double, int): done" << endl;
}

void MCnucl::dumpdNdyTable4Col(char filename[], double *** dNdyTable, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::app);

  double y = rapMin+(rapMax-rapMin)/binRapidity*iy;

    for(int i=0;i<tmax;i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for(int j=0;j<tmax;j++) { // loop over targ thickness
        double ta2 = dT*j;
        if(i>0 && j>0) {
          of << fixed << setprecision(3) << setw(10) <<  y
             << fixed << setprecision(3) << setw(10) <<  ta1
             << fixed << setprecision(3) << setw(10) <<  ta2
             << setprecision(12) << setw(22) <<  dNdyTable[iy][i][j]
             << endl;
        }
      }
    }
     of.close();
}

void MCnucl::dumpdNdydptTable5Col(char filename[], double **** dNdydptTable, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  double y = rapMin+(rapMax-rapMin)/binRapidity*iy;
  double iptmax = MaxPT;

    for(int i=0;i<tmax;i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for(int j=0;j<tmax;j++) { // loop over targ thickness
        double ta2 = dT*j;
        if(i>0 && j>0) {
          for(int ipt=0;ipt<iptmax; ipt++) { //loop over Pt
            double ptstep = PTmin+ipt*dpt;
          of << fixed << setprecision(3) << setw(10) <<  y
             << fixed << setprecision(3) << setw(10) <<  ta1
             << fixed << setprecision(3) << setw(10) <<  ta2
             << fixed << setprecision(3) << setw(10) <<  ptstep;
          of << scientific << setprecision(12) << setw(22) <<  dNdydptTable[iy][i][j][ipt]
             << endl;
             }
        }
      }
    }
     of.close();
}


void MCnucl::deleteNucleus()
{
  proj->clearNucleons();
  targ->clearNucleons();
  
  // Participants are deleted, as they
  // reference the same particles as contained in
  // proj/targ
  
  for(int i=0;i<(int)spectators.size();i++) {
    delete spectators[i];
  }
  for(int i=0;i<(int)binaryCollision.size();i++) {
    delete binaryCollision[i];
  }
  spectators.clear();
  binaryCollision.clear();
}


double MCnucl::Angle(const double x,const double y)
{
    double angl=0.0;
    double r=sqrt(x*x+y*y);
    if(r < 1e-20) return angl;

    if(abs(x)/r < 0.8) {
        //angl=sign(acos(x/r),y)
        angl = y>0 ? abs(acos(x/r)): -abs(acos(x/r));
    }else {
        angl=asin(y/r);
        if(x < 0.0 && angl >= 0.0)
          angl=M_PI-angl;
        else if(x < 0.0)
          angl=-M_PI-angl;
    }

    return angl;
}

void MCnucl::recenterGrid(int iy, int n)
{
  rho->getCMAngle(iy, n);
  rho->recenterPoints(proj->getParticipants(), iy);
  rho->recenterPoints(targ->getParticipants(), iy);
  rho->recenterPoints(binaryCollision, iy);
}

void MCnucl::rotateGrid(int iy, int n)
{
  recenterGrid(iy,n);
  
  rho->rotatePoints(proj->getParticipants(), iy);
  rho->rotatePoints(targ->getParticipants(), iy);
  rho->rotatePoints(binaryCollision, iy);
}


void MCnucl::dumpBinaryTable(char filename[])
{
  double x,y;
  ofstream of;

  of.open(filename, std::ios_base::app);
  for (int idx=0; idx<binaryCollision.size(); idx++)
  {
    x = binaryCollision[idx]->x;
    y = binaryCollision[idx]->y;
    of  << setprecision(3) << setw(10) << x
        << setprecision(3) << setw(10) << y
        << endl;
  }
  of.close();
  
  of.open("data/wounded.data");
  proj->dumpParticipants(of);
  targ->dumpParticipants(of);
  of.close();

  of.open("data/nucl1.data");
  proj->dumpNucleons(of);
  of.close();

  of.open("data/nucl2.data");
  targ->dumpNucleons(of);
  of.close();
  
}

void MCnucl::dumpparticipantTable(char filename[])
{
  ofstream of;
  of.open(filename, std::ios_base::app);
  proj->dumpParticipants(of);
  targ->dumpParticipants(of);
  of.close();
}


int MCnucl::getSpectators()
{
  int count = 0;
  double ecm = paraRdr->getVal("ecm");
  //calculate the rapidity_Y for spectators at a given collision energy
  double v_z = sqrt(1. - 1./((ecm/2.)*(ecm/2.)));
  double rapidity_Y = 0.5*log((1. + v_z)/(1. - v_z + 1e-100));
  vector<Particle*> nucl1 = proj->getNucleons();
  vector<Particle*> nucl2 = targ->getNucleons();
  for(int i=0;i<(int)nucl1.size();i++) { // loop over proj. nucleons
    if(nucl1[i]->getNumberOfCollision()==0) {
      double x1 = nucl1[i]->x;
      double y1 = nucl1[i]->y;
      spectators.push_back(new Spectator(x1, y1, rapidity_Y));
      count++;
    }
  }
  for(int i=0;i<(int)nucl2.size();i++) { // loop over targ. nucleons
    if(nucl2[i]->getNumberOfCollision()==0) {
      double x2 = nucl2[i]->x;
      double y2 = nucl2[i]->y;
      spectators.push_back(new Spectator(x2, y2, -rapidity_Y));
      count++;
    }
  }
  return(count);
}

void MCnucl::dumpSpectatorsTable(int event)
{
  double x, y, rap;
  ostringstream of_stream; 
  ofstream of;
  of_stream << "data/Spectators_event_" << event << ".dat";
  
  of.open(of_stream.str().c_str());
  for (int idx=0; idx<spectators.size(); idx++)
  {
    x = spectators[idx]->getX();
    y = spectators[idx]->getY();
    rap = spectators[idx]->getRapidity_Y();
    of  << scientific << setprecision(4) << setw(10) 
        << x << "  " << y << "  " << rap
        << endl;
  }
  of.close();
}

double MCnucl::sampleFluctuationFactorforParticipant()
{
   double eps = 1e-8;
   double fluctfactor = 1.0;
   double Gamma_k = 1./ccFluctuationGammaTheta;
   double k_part = (1 - Alpha + eps)/2.*Gamma_k;
   double theta_part = 2./(1 - Alpha + eps)*ccFluctuationGammaTheta;
   if(CCFluctuationModel == 6)  //Gamma distribution for MC-Glauber
      fluctfactor = gsl_ran_gamma(gslRng, k_part, theta_part);
   
   return(fluctfactor);
}

double MCnucl::sampleFluctuationFactorforBinaryCollision()
{
   double eps = 1e-8;
   double fluctfactor = 1.0;
   double Gamma_k = 1./ccFluctuationGammaTheta;
   double k_binary = (Alpha+eps)*Gamma_k;
   double theta_binary = 1./(Alpha+eps)*ccFluctuationGammaTheta;
   if(CCFluctuationModel == 6)  //Gamma distribution for MC-Glauber
      fluctfactor = gsl_ran_gamma(gslRng, k_binary, theta_binary);
   
   return(fluctfactor);
}
