// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <dirent.h>
#include <ctime>
#include <bitset>
#include <cmath>
#include <chrono>
#include <time.h>
#include <vector>
#include <queue>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <queue>
#include <iterator>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TGraph.h>
#include <GeoOctuplet.hh>
#include <Hit.hh>
#include <Road.hh>
#include <TStyle.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <unistd.h>
#include <time.h>
#include "VectorDict.cxx"

using namespace std;

TRandom3 *ran = new TRandom3(time(NULL));

bool db = true; // debug output flag

// SOME CONSTANTS

int NPLANES = 8;
int NSTRIPS;
int ADDC_BUFFER = 8;
double xlow, xhigh, ylow, yhigh;             // chamber dimensions
double mu_xlow, mu_xhigh, mu_ylow, mu_yhigh; // active chamber area to decouple edge effects

int XROAD, UVFACTOR;

int NSTRIPS_UP_UV, NSTRIPS_DN_UV;
int NSTRIPS_UP_XX, NSTRIPS_DN_XX;

int bc_wind;
int sig_art;

double B = (1 / TMath::Tan(1.5 / 180. * TMath::Pi()));

bool compare_age(Hit *a, Hit *b)
{
   return (a->Age() < b->Age());
}

bool compare_channel(Hit *a, Hit *b)
{
   return (a->Channel() < b->Channel());
}

bool compare_second(std::pair<int, double> a, std::pair<int, double> b)
{
   return (a.second < b.second);
}

void set_chamber(string chamber, int m_wind, int m_sig_art, int m_xroad, bool uvrflag, bool trapflag, int m_nstrips)
{
   // function to set parameters in a smart way

   if (chamber == "small") {
      if (m_nstrips == -1) {
         NSTRIPS = 8800; // has to be multiple of x road
      } else {
         NSTRIPS = m_nstrips;
      }
      xlow = 0.;
      xhigh = NSTRIPS * 0.4 - 0.2;
      ylow = 0.;
      yhigh = trapflag ? 1821. : 500.;
   } else if (chamber == "large") {
      if (m_nstrips == -1) {
         NSTRIPS = 8800; // has to be multiple of x road
      } else {
         NSTRIPS = m_nstrips;
      }
      xlow = 0.;
      xhigh = NSTRIPS * 0.4 - 0.2;
      ylow = 0.;
      yhigh = 2200.;
   } else if (chamber == "oct") {
      if (m_nstrips == -1) {
         NSTRIPS = 512; // has to be multiple of x road
      } else {
         NSTRIPS = m_nstrips;
      }
      xlow = 0.;
      xhigh = NSTRIPS * 0.4 - 0.2;
      ylow = 17.9;
      yhigh = 217.9;
   } else {
      cout << "Now i Am here " << endl;
      exit(EXIT_FAILURE);
   }

   // active area
   mu_xlow = 100 * 0.4 + 0.2;
   mu_xhigh = NSTRIPS * 0.4 - 0.2 - 100 * 0.4;

   mu_ylow = ylow;
   mu_yhigh = yhigh;

   bc_wind = m_wind;

   sig_art = m_sig_art;

   XROAD = m_xroad;
   UVFACTOR = round((yhigh - ylow) / (B * 0.4 * 2) / XROAD);

   if (!uvrflag) {
      // this is for 8 strip x-roads, i think
      NSTRIPS_UP_UV = UVFACTOR * XROAD + NSTRIPS_UP_XX;
      NSTRIPS_DN_UV = UVFACTOR * XROAD;
      NSTRIPS_UP_XX = 4;
      NSTRIPS_DN_XX = 0;

   } else {
      NSTRIPS_UP_UV = 4;
      NSTRIPS_DN_UV = 0;
      NSTRIPS_UP_XX = 4;
      NSTRIPS_DN_XX = 0;
   }
}

////<--------------CHECKS IF PARAMS IN RANGE--------------->

int inside_trapezoid(double x, double y, double inner_radius, double outer_radius, double length, double base_width,
                     double top_width)
{
   double slope = (outer_radius - inner_radius) / ((top_width - base_width) / 2.0);
   double offset = inner_radius - (slope * base_width / 2.0);
   if (x > outer_radius)
      return 0; // top
   if (x < inner_radius)
      return 0; // bottom
   if (y > (x - offset) / slope)
      return 0; // right
   if (y < -(x - offset) / slope)
      return 0; // left
   return 1;
}

int fiducial(double x, double y, string chamber)
{

   if (chamber != "large" && chamber != "small") {
      cerr << "fiducial doesnt understand this chamber: " << chamber << endl;
      return -1;
   }
   int large = (chamber == "large");

   double NSW_MM1_InnerRadius = 0; // chambers starts at 0 // large ? 923.0  : 895.0;
   double NSW_MM1_Length = large ? 2310.0 : 2210.0;
   double NSW_MM2_Length = large ? 1410.0 : 1350.0;
   double NSW_MM1_baseWidth = large ? 640.0 : 500.0;
   double NSW_MM1_topWidth = large ? 2008.5 : 1319.2;
   double NSW_MM2_baseWidth = large ? 2022.8 : 1321.1;
   double NSW_MM2_topWidth = large ? 2220.0 : 1821.5;
   double NSW_MM1_outerRadius = NSW_MM1_InnerRadius + NSW_MM1_Length;
   double NSW_MM2_InnerRadius = NSW_MM1_outerRadius;
   double NSW_MM2_outerRadius = NSW_MM2_InnerRadius + NSW_MM2_Length;

   if (inside_trapezoid(x, y - (yhigh + ylow) / 2.0, NSW_MM1_InnerRadius, NSW_MM1_outerRadius, NSW_MM1_Length,
                        NSW_MM1_baseWidth, NSW_MM1_topWidth))
      return 1;
   if (inside_trapezoid(x, y - (yhigh + ylow) / 2.0, NSW_MM2_InnerRadius, NSW_MM2_outerRadius, NSW_MM2_Length,
                        NSW_MM2_baseWidth, NSW_MM2_topWidth))
      return 1;
   return 0;
}

///<--------RETURNS THX, THY ---------------->
tuple<double, double> cosmic_angle(int angcos, double angx, double angy)
{
   double thx = ran->Uniform(-angx, angx) * TMath::Pi() / 180;
   double thy = ran->Uniform(-angy, angy) * TMath::Pi() / 180;
   return make_tuple(thx, thy);
}

// <<<---------Calculates RATE 'rate'------------------>
double predicted_rate(int strip, string chamber)
{

   if (chamber != "large" && chamber != "small") {
      cerr << "predicted_rate doesnt understand this chamber: " << chamber << endl;
      return -1.0;
   }
   int large = (chamber == "large");
   double pitch = 0.4;
   double offset = large ? 923.0 : 895.0;
   double r = offset + pitch * (double)(strip);

   // mm->cm
   r = r / 10;

   double rate = 0.0;
   if (large)
      rate = (-9.938824) + (6288.351422) / r + (45942.902843) / pow(r, 2);
   else
      rate = (-5.018321) + (3396.877744) / r + (164524.202988) / pow(r, 2);

   // kHz->Hz
   return rate * 1000;
}

//<---------------_GENERATES MUONS_ tuple (x, y, thx, thy)---------->
tuple<double, double, double, double> generate_muon(vector<double> &xpos, vector<double> &ypos, vector<double> &zpos,
                                                    string chamber, int angcos, double angx, double angy, bool trapflag)
{

   double x = 9e9;
   double y = 9e9;

   if (trapflag) {
      while (!fiducial(x, y, chamber)) {
         x = ran->Uniform(mu_xlow, mu_xhigh);
         y = ran->Uniform(mu_ylow, mu_yhigh);
      }
   } else {
      x = ran->Uniform(mu_xlow, mu_xhigh);
      y = ran->Uniform(mu_ylow, mu_yhigh);
   }

   double thx, thy;

   std::tie(thx, thy) = cosmic_angle(angcos, angx, angy);

   double avgz = 0.5 * (zpos[0] + zpos[NPLANES - 1]);
   double x_b, y_b;
   // double z, x_b, y_b;
   for (int j = 0; j < NPLANES; j++) {
      // z = zpos[j];
      x_b = TMath::Tan(thx) * (zpos[j] - avgz) + x;
      y_b = TMath::Tan(thy) * (zpos[j] - avgz) + y;
      xpos[j] = x_b;
      ypos[j] = y_b;
   }
   return make_tuple(x, y, thx, thy);
}

//<-----------------------GENERATES BACKGROUNG (vector 'bkghits')-----------------______>>>>>>>>>

vector<Hit *> generate_bkg(int start_bc, const GeoOctuplet &geometry, int bkgrate, string chamber)
{

   vector<Hit *> bkghits;

   int noise_window = bc_wind * 5;
   int start_noise = start_bc - bc_wind * 2; // this takes into account overwriting tp hits for a plane+road

   // assume uniform distribution of background - correct for noise
   double time_window = noise_window * 25 * pow(10, -9);
   double bkg_prob = bkgrate * time_window;
   for (int j = 0; j < NPLANES; j++) {
      // int nbkg = expbkg;
      for (int k = 0; k < NSTRIPS; k++) {
         double prob = ran->Uniform(0, 1.);
         if (bkgrate == -1)
            bkg_prob = predicted_rate(k, chamber) * time_window;
         if (prob < bkg_prob) {
            Hit *newhit = nullptr;
            newhit = new Hit(j, start_noise + ran->Integer(noise_window), k, true, geometry);
            bkghits.push_back(newhit);
         }
      }
   }
   return bkghits;
}

//<---CREATES A VECTOR OF PLANES WHICH REGISTERED HITS    ('oct_hitmask')------->>>>>
vector<int> oct_response(vector<double> &xpos, vector<double> &ypos, vector<double> &zpos, vector<double> &mm_eff)
{
   // gives detector response to muon, returns list of which planes registered hit

   int n_mm = 0;
   vector<int> oct_hitmask(NPLANES, 0);
   for (int j = 0; j < NPLANES; j++) {
      if (ran->Uniform(0., 1.) < mm_eff[j]) {
         oct_hitmask[j] = 1;
         n_mm++;
      }
   }
   return oct_hitmask;
}

// Prints Binary for ANY datatype
// assumes little endian
void printBits(size_t const size, void const *const ptr)
{
   unsigned char *b = (unsigned char *)ptr;
   unsigned char byte;
   int i, j;

   for (i = size - 1; i >= 0; i--) {
      for (j = 7; j >= 0; j--) {
         byte = (b[i] >> j) & 1;
         printf("%u", byte);
      }
   }
   puts("");
}

#pragma pack(1)
// union Pack {
//    unsigned int code;
//    struct {
//       unsigned int vmm : 4;
//       unsigned int channel : 7; // 5 to reset it 0-63
//       unsigned int addc : 4;
//       unsigned int mmfe : 3;
//       unsigned int age : 5; // reset when defining BC
//    } storage_t;
// };
//<-------------TO DO TO DO TO DOOOOOOOOOOOOOOOOOOOOOOOOOOO------------------>
struct storage_t {
   unsigned int vmm; //
   unsigned int channel;
   unsigned int addc;
   unsigned int mmfe;
   unsigned int age; // total of 23 bits , prints 24
};

// unsigned int pack(unsigned int v, unsigned int c, unsigned int a, unsigned int m, unsigned int ag)
// {
//    unsigned int code = c;
//    code *= ;
//    code += c;
//    code *= 100;
//    code += b;
//    code *= 18;
//    code += a;
//    return code;
// }

vector<storage_t> finder(vector<Hit *> hits, int mu_firstbc, bool ideal_vmm, bool ideal_addc, bool ideal_tp, int nevts)
{
   // applies the MMTP finder to a series of hits and roads
   // returns slope structs for roads which found a coincidence and have at least 1 real muon hit
   int bc_start = 999999;
   int bc_end = -1;
   unsigned int ibc = 0;
   vector<storage_t> storage = {};

   if (hits.size() == 0) {
      cout << "Now i Am here " << endl;
      return storage;
   }

   std::sort(hits.begin(), hits.end(), compare_age);
   bc_start = hits.front()->Age();
   bc_end = hits.back()->Age();
   bc_start = bc_start - bc_wind * 2;
   bc_end = bc_end + bc_wind * 2;

   vector<Hit *> hits_now = {};
   vector<int> vmm_same = {};
   vector<pair<int, double>> addc_same = {};
   vector<int> to_erase = {};
   int n_vmm = NSTRIPS / 64;
   int n_addc = NSTRIPS / 2048;

   // each road makes independent triggers, evaluated on each BC
   for (int bc = bc_start; bc < bc_end; bc++) {

      hits_now.clear();

      ///!!!!!!!!!!!!!!!!!< ---------Asssigns hits by BC ---------------------->>>>.

      for (unsigned int j = ibc; j < hits.size(); j++) { ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! very important function
         if (hits[j]->Age() == bc) {
            hits_now.push_back(hits[j]);
         } else if (hits[j]->Age() > bc) {
            ibc = j;
            break;
         }
      }

      ///!!!!!!!!!!!!!!!!!< ---------Groups hits by MMFE8Index and keeps 1 hit per chip---------------------->>>>.

      // implement vmm ART-like signal
      for (int ib = 0; ib < NPLANES; ib++) {
         for (int j = 0; j < n_vmm; j++) {
            vmm_same.clear();
            // save indices of all elements in vmm j
            for (unsigned int k = 0; k < hits_now.size(); k++) {
               if (hits_now[k]->MMFE8Index() != ib)
                  continue;
               if (hits_now[k]->VMM() == j) {
                  vmm_same.push_back(k);
               }
            }
            // if 2+ hits in same vmm, erase all except 1 randomly
            if (vmm_same.size() > 1) {
               int the_chosen_one = ran->Integer((int)(vmm_same.size()));

               if (ideal_vmm)
                  for (unsigned int k = 0; k < vmm_same.size(); k++)
                     if (hits_now[vmm_same[k]]->IsReal()) {
                        the_chosen_one = k;
                        break;
                     }
               for (int k = vmm_same.size() - 1; k > -1; k--)
                  if (k != the_chosen_one)
                     hits_now.erase(hits_now.begin() + vmm_same[k]);
            }
         }

         // implement ADDC-like handling
         for (int ia = 0; ia < n_addc; ia++) {

            // collect all hits on one ADDC
            addc_same.clear();
            for (unsigned int k = 0; k < hits_now.size(); k++)
               if (hits_now[k]->MMFE8Index() == ib && hits_now[k]->ADDC() == ia)
                  addc_same.push_back(std::make_pair(k, hits_now[k]->Channel()));

            if ((int)(addc_same.size()) <= ADDC_BUFFER)
               continue;

            // priority encode the hits by channel number; remember hits 8+
            to_erase.clear();
            std::sort(addc_same.begin(), addc_same.end(), compare_second);
            for (int it = ADDC_BUFFER; it < (int)(addc_same.size()); it++)
               to_erase.push_back(addc_same[it].first);

            // reverse and erase
            std::sort(to_erase.rbegin(), to_erase.rend());
            for (auto k : to_erase)
               if (ideal_addc && hits_now[k]->IsReal())
                  continue;
               else
                  hits_now.erase(hits_now.begin() + k);
         }
      }
      for (unsigned int m = 0; m < hits_now.size(); m++) {
         storage_t m_storage;
         m_storage.vmm = round((hits_now[m]->VMM()) / 18);
         m_storage.channel = round((hits_now[m]->Channel()) / 138);
         m_storage.mmfe = hits_now[m]->MMFE8Index();
         m_storage.addc = hits_now[m]->ADDC();
         //    m_storage.age = hits_now[m]->Age();
         m_storage.age = ran->Uniform(0, nevts);
         storage.push_back(m_storage);
      }
   }
   return storage;
}

//    /// FIFO functions
//    class Queue {
//    private:
//       int A[MAX_SIZE];
//       int front, rear;

//    public:
//       // Constructor - set front and rear as -1.
//       // We are assuming that for an empty Queue, both front and rear will be -1.
//       Queue()
//       {
//          front = -1;
//          rear = -1;
//       }

//       // To check wheter Queue is empty or not
//       bool IsEmpty() { return (front == -1 && rear == -1); }

//       // To check whether Queue is full or not
//       bool IsFull() { return (rear + 1) % MAX_SIZE == front ? true : false; }

//       // Inserts an element in queue at rear end
//       void Enqueue(int x)
//       {
//          cout << "Enqueuing " << x << " \n";
//          if (IsFull()) {
//             cout << "Error: Queue is Full\n";
//             return;
//          }
//          if (IsEmpty()) {
//             front = rear = 0;
//          } else {
//             rear = (rear + 1) % MAX_SIZE;
//          }
//          A[rear] = x;
//       }

//       // Removes an element in Queue from front end.
//       void Dequeue()
//       {
//          cout << "Dequeuing \n";
//          if (IsEmpty()) {
//             cout << "Error: Queue is Empty\n";
//             return;
//          } else if (front == rear) {
//             rear = front = -1;
//          } else {
//             front = (front + 1) % MAX_SIZE;
//          }
//       }
//       // Returns element at front of queue.
//       int Front()
//       {
//          if (front == -1) {
//             cout << "Error: cannot return front from empty queue\n";
//             return -1;
//          }
//          return A[front];
//       }
//       /*
//          Printing the elements in queue from front to rear.
//          This function is only to test the code.
//          This is not a standard function for Queue implementation.
//       */
//       void Print()
//       {
//          // Finding number of elements in queue
//          int count = (rear + MAX_SIZE - front) % MAX_SIZE + 1;
//          cout << "Queue       : ";
//          for (int i = 0; i < count; i++) {
//             int index = (front + i) % MAX_SIZE; // Index of element while travesing circularly from front
//             cout << A[index] << " ";
//          }
//          cout << "\n\n";
//       }
//    };

unsigned int
pack(unsigned int p_vmm, unsigned int p_channel, unsigned int p_addc, unsigned int p_mmfe, unsigned int p_age)
{
   unsigned int code = p_age;
   code *= 64;
   code += p_channel;
   code *= 8;
   code += p_vmm;
   code *= 8;
   code += p_mmfe;
   code *= 5;
   code += p_addc;
   return code;
}
unsigned int unpack(unsigned int code, unsigned int p_vmm, unsigned int p_mmfe, unsigned int p_channel,
                    unsigned int p_addc, unsigned int p_age)
{
   p_addc = code % 5;
   code /= 5;
   p_mmfe = code % 8;
   code /= 8;
   p_vmm = code % 8;
   code /= 8;
   p_channel = code % 64;
   code /= 64;
   p_age = code;
   return p_age;
}
using namespace std;

void showq(queue<int> gq)
{
   queue<int> g = gq;
   while (!g.empty()) {
      cout << '\t' << g.front();
      g.pop();
   }
   cout << '\n';
}
double memory = 0;

int main(int argc, char *argv[])
{
   int nevents = -1;
   int bkgrate = 0; // Hz per strip

   // #define MAX_SIZE nevents

   //       Queue Q; // creating an instance of Queue.
   //       Q.Enqueue(2);
   //       Q.Print();
   //       Q.Enqueue(4);
   //       Q.Print();
   //       Q.Enqueue(6);
   //       Q.Print();
   //       Q.Dequeue();
   //       Q.Print();
   //       Q.Enqueue(8);
   //       Q.Print();

   int m_xroad = 8;
   int m_NSTRIPS = -1;

   int m_bcwind = 8;
   int m_sig_art = 32;
   int m_sig_art_x = 32;
   vector<double> mm_eff = {1., 1., 1., 1., 1., 1., 1., 1.};
   double chamber_eff = -1;

   double angx = 0;
   double angy = 0;
   int angcos = 0;
   // coincidence params
   int m_xthr = 2;
   int m_uvthr = 2;

   bool bkgflag = false;
   bool uvrflag = false;
   bool trapflag = false;
   bool ideal_tp = false;
   bool ideal_vmm = false;
   bool ideal_addc = false;
   bool smear_art = false;
   bool funcsmear_art = false;
   bool bkgonly = false;
   bool write_tree = false;

   char chamberType[400];
   char outputFileName[400];

   if (argc < 3) {
      cout << "Error at Input: please specify number of events to generate " << endl;
      cout << "Example:   ./sim -n 100 -ch <chamber type> -o output.root" << endl;
      cout << "Example:   ./sim -n 100 -ch <chamber type> -b <bkg rate in kHz/strip> -o output.root" << endl;
      cout << "Example:   ./sim -n 100 -ch <chamber type> -b <bkg rate in Hz/strip> -p <make event displays> -o "
              "output.root"
           << endl;
      cout << "Other options include: -w <bc_wind> -sig <art res (ns)>" << endl;
      cout << "If art res = 0, then we do bkg only" << endl;
      return 0;
   }

   bool b_out = false;
   bool ch_type = false;
   for (int i = 1; i < argc; i++) {
      if (strncmp(argv[i], "-n", 2) == 0) {
         nevents = atoi(argv[i + 1]); /// nevents
      }
      if (strncmp(argv[i], "-o", 2) == 0) {
         sscanf(argv[i + 1], "%s", outputFileName);
         b_out = true; // bout
      }
      if (strncmp(argv[i], "-x", 2) == 0) {
         m_xroad = atoi(argv[i + 1]); // m_xroad
      }
      if (strncmp(argv[i], "-w", 2) == 0) {
         m_bcwind = atoi(argv[i + 1]); // m_bcwind
      }
      if (strncmp(argv[i], "-thrx", 5) == 0) { // m
         m_xthr = atoi(argv[i + 1]);
      }
      if (strncmp(argv[i], "-thruv", 6) == 0) {
         m_uvthr = atoi(argv[i + 1]);
      }
      if (strncmp(argv[i], "--trap", 6) == 0) {
         trapflag = true;
      }
      if (strncmp(argv[i], "-ch", 3) == 0) {
         sscanf(argv[i + 1], "%s", chamberType);
         ch_type = true;
      }
      if (strncmp(argv[i], "-b", 2) == 0) {
         bkgrate = atoi(argv[i + 1]);
         bkgflag = true;
      }
      if (strncmp(argv[i], "-e", 2) == 0) {
         chamber_eff = atof(argv[i + 1]);
         for (unsigned int i = 0; i < mm_eff.size(); i++)
            mm_eff[i] = chamber_eff;
      }
      if (strncmp(argv[i], "-angx", 5) == 0) {
         angx = fabs(atof(argv[i + 1]));
      }
      if (strncmp(argv[i], "-angy", 5) == 0) {
         angy = fabs(atof(argv[i + 1]));
      }
      if (strncmp(argv[i], "-angcos", 7) == 0) {
         angcos = 1;
      }
      if (strncmp(argv[i], "-tree", 5) == 0) {
         write_tree = true;
      }
      if (strncmp(argv[i], "-ideal-vmm", 10) == 0) {
         ideal_vmm = true;
      }
      if (strncmp(argv[i], "-ideal-addc", 11) == 0) {
         ideal_addc = true;
      }
      if (strncmp(argv[i], "-ideal-tp", 9) == 0) {
         ideal_tp = true;
      }
      if (strncmp(argv[i], "-strips", 7) == 0) {
         m_NSTRIPS = atoi(argv[i + 1]);
      }
      if (strncmp(argv[i], "-smearstrips", 12) == 0) {
         m_sig_art_x = atoi(argv[i + 1]);
         smear_art = true;
      }
      if (strncmp(argv[i], "-smear", 6) == 0) {
         smear_art = true;
      }
      if (strncmp(argv[i], "-funcsmear", 10) == 0) {
         funcsmear_art = true;
      }
   }

   if (!b_out) {
      cout << "Error at Input: please specify output file (-o flag)" << endl;
      return 0;
   }

   if (!ch_type) {
      cout << "Error at Input: please specify chamber type (-ch flag, options: large, small, oct)" << endl;
      return 0;
   }

   if (nevents == -1) {
      cout << "Didn't set the number of generated events! Exiting." << endl;
      return 0;
   }
   set_chamber(string(chamberType), m_bcwind, m_sig_art, m_xroad, uvrflag, trapflag, m_NSTRIPS);
   //    void set_chamber(string chamber, int m_wind, int m_sig_art, int m_xroad, bool uvrflag, bool trapflag, int
   //    m_nstrips)

   if (NSTRIPS % XROAD != 0) {
      cout << "Number of strips not divisible by the road size!" << endl;
      return 0;
   }

   if (((mu_xlow || mu_xhigh) < xlow) || ((mu_xlow || mu_xhigh) > xhigh) || (mu_xlow > mu_xhigh)) {
      cout << "Muon active area is outside the chamber area!" << endl;
      return 0;
   }

   cout << "--------------" << endl;
   cout << "OCT SIM ✪ ‿ ✪ " << endl;
   cout << "--------------" << endl;

   printf("\r >> x-road size (in strips): %d, +/- neighbor roads (uv): %d", XROAD, UVFACTOR);
   cout << endl;
   printf("\r >> art res (in ns): %d", m_sig_art);
   cout << endl;
   cout << "\r >> Using BCID window: " << bc_wind << endl;
   printf("\r >> Background rate of %d Hz per strip", bkgrate);
   cout << endl;
   printf("\r >> Assuming chamber size: (%4.1f,%4.1f) in mm", xhigh - xlow, yhigh - ylow);
   cout << endl;
   printf("\r >> Assuming muon active area: (%4.1f,%4.1f) in mm", mu_xhigh - mu_xlow, mu_yhigh - mu_ylow);
   cout << endl;
   printf("\r >> Using UV roads: %s", (uvrflag) ? "true" : "false");
   cout << endl;
   printf("\r >> Using trapezoidal geometry: %s", (trapflag) ? "true" : "false");
   cout << endl;
   printf("\r >> Using thresholds (x, uv): (%d, %d)", m_xthr, m_uvthr);
   cout << endl;
   printf("\r >> Generate muons with cosmic distribution: %s", (angcos) ? "true" : "false");
   cout << endl;
   printf("\r >> Generate muons with angle (x) from %f to %f", -angx, angx);
   cout << endl;
   printf("\r >> Generate muons with angle (y) from %f to %f", -angy, angy);
   cout << endl;

   cout << "Generating " << nevents << " events" << endl;

   // define output ntuple

   vector<vector<int>> Hit_strips;
   vector<vector<int>> Hit_planes;
   vector<vector<int>> Hit_ages;
   vector<int> N_muon;
   vector<int> N_xmuon;

   vector<double> trig_x;
   vector<double> trig_y;

   // geometry stuff
   double xlen = xhigh - xlow;
   double ylen = yhigh - ylow;

   GeoOctuplet *GEOMETRY;
   if (string(chamberType) == "oct")
      GEOMETRY = new GeoOctuplet(false, xlen, ylen);
   else
      GEOMETRY = new GeoOctuplet(true, xlen, ylen);

   //    queue<int> fifo;

   for (int i = 0; i < nevents; i++) {

      // generate muon
      double co = 2.7;
      co = 0;
      vector<double> zpos = {-co, 11.2 + co, 32.4 - co, 43.6 + co, 113.6 - co, 124.8 + co, 146.0 - co, 157.2 + co};
      vector<double> xpos(NPLANES, -1.);
      vector<double> ypos(NPLANES, -1.);

      double xmuon, ymuon, thx, thy;
      std::tie(xmuon, ymuon, thx, thy) =
         generate_muon(xpos, ypos, zpos, string(chamberType), angcos, angx, angy, trapflag);

      double real_x_muon;
      double real_y_muon;
      real_x_muon = xmuon;
      real_y_muon = ymuon;

      vector<int> oct_hitmask = oct_response(xpos, ypos, zpos, mm_eff);
      vector<int> art_bc(NPLANES, -1.);
      double smallest_bc = 999999.;

      vector<Hit *> hits;
      cout << "This is the size of hits originally " << hits.size() << endl;

      int n_u = 0;
      int n_v = 0;
      int n_x1 = 0;
      int n_x2 = 0;

      double art_time;

      double strip, strip_smear;

      for (int j = 0; j < NPLANES; j++) {
         if (oct_hitmask[j] == 1) {
            if (j < 2)
               n_x1++;
            else if (j > 5)
               n_x2++;
            else if (j == 2 || j == 4)
               n_u++;
            else
               n_v++;
            art_time = ran->Gaus(400., (double)(sig_art));
            art_bc[j] = (int)floor(art_time / 25.);
            Hit *newhit = nullptr;

            strip = GEOMETRY->Get(j).channel_from_pos(xpos[j], ypos[j]);
            if (smear_art) {
               strip_smear = round(ran->Gaus(strip, m_sig_art_x));
            } else {
               strip_smear = strip;
            }
            newhit = new Hit(j, art_bc[j], strip_smear, false, *GEOMETRY);
            // newhit = new Hit(j, art_bc[j], xpos[j], ypos[j], false, *GEOMETRY);
            if (!bkgonly)
               hits.push_back(newhit);
         }
      }
      cout << "HITS after oct_hitmask before finder:hits[1] " << hits[1] << endl;

      if (true) {
         cout << "N muonhits: " << hits.size() << endl;
         for (unsigned int j = 0; j < hits.size(); j++) {
            printf("Muon hit (mmfe8, BC, strip, addc): (%d,%d,%4.1f,%d)\n", hits[j]->MMFE8Index(), hits[j]->Age(),
                   hits[j]->Channel(), hits[j]->ADDC());
         }
      }

      for (unsigned int j = 0; j < art_bc.size(); j++) {
         if (art_bc[j] == -1)
            continue;
         else if (art_bc[j] < smallest_bc)
            smallest_bc = art_bc[j];
      }

      // assume bkg rate has oct_response factored in!

      vector<Hit *> all_hits = hits;

      if (bkgflag) {
         vector<Hit *> bkghits = generate_bkg(smallest_bc, *GEOMETRY, bkgrate, string(chamberType));
         if (db)
            cout << "Number of bkg hits: " << bkghits.size() << endl;
         all_hits.insert(all_hits.end(), bkghits.begin(), bkghits.end());
      }

      //    for (unsigned int ihit = 0; ihit < all_hits.size(); ihit++) {
      //       int ib = all_hits[ihit]->MMFE8Index();
      //       printf("%d , ", ib);
      //       //    if (all_hits[ihit]->IsNoise() && (ib < 2 || ib > 5))
      //       //       printf("%f\n",
      //       //              GEOMETRY->Get(ib).LocalXatYend(all_hits[ihit]->Channel()) +
      //       GEOMETRY->Get(ib).Origin().X() -
      //       //              xmuon);
      //    }

      // if (db)
      //    cout << "Total number of muons hit+background " << all_hits.size() << endl;
      vector<storage_t> outstorage;
      outstorage = finder(all_hits, smallest_bc, ideal_vmm, ideal_addc, ideal_tp, nevents);
      // cout << "Total number of hits entirely " << outstorage.size() << endl;

      cout << endl;
      cout << "--------------" << endl;
      cout << "OCT SIM ✪ ‿ ✪ " << endl;
      cout << "--------------" << endl;
      for (unsigned int j = 0; j < outstorage.size(); j++) {
         // declaring output string stream
         //    ostringstream str1;

         //    // Sending a number as a stream into output
         //    // string
         //    str1 << outstorage[j].addc;
         //    str1 << outstorage[j].mmfe;
         //    str1 << round((outstorage[j].vmm) / 17.1875);
         //    str1 << round((outstorage[j].channel) / 137.5);
         //    str1 << outstorage[j].age;

         //    // the str() coverts number into string
         //    string mystring = str1.str();
         //    cout << "STRING :  " << mystring << endl;
         //    cout << " ADDC " << outstorage[j].addc << " MMFE8Index " << outstorage[j].mmfe << "  vmm  "
         //         << (outstorage[j].vmm) << "  channel " << (outstorage[j].channel) << " Age: " << outstorage[j].age
         //         << endl;
         unsigned int myoutput =
            pack(outstorage[j].vmm, outstorage[j].channel, outstorage[j].addc, outstorage[j].mmfe, outstorage[j].age);
         memory = memory + sizeof(myoutput);
         //    fifo.push(myoutput);
         //    cout << "SIZE OF  my string:  " << sizeof(mystring) << endl;
         //    str1.str(std::string());
      }

      // unsigned int newvmm;
      // newvmm = round((outstorage[0].vmm) / 17.1875);

      // cout << "BINARY REPRESENTATION simple vmm  " << endl;

      // printBits(sizeof(outstorage[0]), &(outstorage[0]));
      // cout << "BINARY REPRESENTATION shifted vmm  " << sizeof(outstorage[0]) << endl;
      // unsigned int myoutput =
      //    pack(outstorage[0].vmm, outstorage[0].channel, outstorage[0].addc, outstorage[0].mmfe, outstorage[0].age);
      // cout << " CODE CODE :  " << myoutput << endl;

      // cout << " CODE CODE :  " << myoutput << endl;
      // cout << " SIZE CODE :  " << sizeof(myoutput) << endl;

      // BINARY REPRESENTATION
      // cout << "  BINARY REPRESENTATION of CODE  " << endl;
      // printBits(sizeof(myoutput), &(myoutput));
      // cout << " BINARY REPRESENTATION of string CODE  " << endl;
      // string mystring = 'std::to_string(myoutput)';
      // printBits(sizeof(mystring), &(mystring));
      // cout << " Size of string CODE  " << sizeof(mystring) << endl;

      // UNPACKING
      // int chunk = 3;
      // unsigned int *buffer;
      // buffer = (unsigned int *)malloc(chunk);
      // buffer = &myoutput;
      // printBits(sizeof(buffer), &(buffer));
      // cout << " BUFFER :  " << &buffer << endl;
      // cout << " Size of BUFFER :  " << sizeof(buffer) << endl;
      // unsigned int testchannel = unpack(myoutput, outstorage[0].vmm, outstorage[0].mmfe, outstorage[0].channel,
      //                                   outstorage[0].addc, outstorage[0].age);
      // cout << " channel from unpack CODE :  " << testchannel << endl;

      // std::string structToString(outstorage[0]);
      // {
      //    char *tempCstring = new char[100];
      //    memcpy(tempCstring, theStruct, 100);
      //    tempCstring[101] = '0';
      //    std::string returnVal(tempCstring, 100);
      //    delete tempCstring;
      //    cout << "STRING:  " << returnVal << endl;
      // }
   }
   vector<Double_t> sizeandbkg = {};

   TFile *fout = new TFile(outputFileName, "RECREATE");
   TTree *tree = new TTree("gingko", "gingko");
   tree->Branch("sizeandbkg", &sizeandbkg);
   sizeandbkg.push_back(memory);
   cout << "memory after n events:  " << memory << endl;

   tree->Fill();
   fout->cd();
   tree->Write();

   fout->cd();
   fout->Close();
}