/* 
 * Kepler Q1-Q17 DR24 Robovetter
 * 
 * Compile via: g++ -std=c++11 -O3 -o robovet DR24-RoboVetter.cpp
 * 
 * Run as ./robovet INFILE  OUTFILE
 * 
 * For example:
 * 
 * "./robovet RoboVetter-Input.txt RoboVetter-Output.txt"   for the real data
 * 
 * "./robovet RoboVetter-Inject-Input.txt RoboVetter-Inject-Output.txt"   for the artifically injected transit data
 *
 */


#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <math.h>

using namespace std;

// Declare Functions
void READDATA(),TRANSITLIKE(),SECECLIPSE(),ISSEC(),COMPPT(double,double,double,double);
double INVERFC(double), PSIG(double,double), ESIG(double,double,double);

// Declare constants
const int NMAX = 50000; // Max number of TCEs. Increase as needed.
const double PSIG_THRESH = 3.25;   // Period matching threshold
const double ESIG_THRESH = 2.0;    // Epoch matching threshold
const double WIDTHFAC = 2.5;       // Transit exclusion width actor
const double MISSIONDUR = 1600.0;  // Mission duration in days

// Declare variables
int i,j,k,l;  // Counting integers
int ntces;  // Number of TCEs to robo-vet
string infilename,outfilename;
ifstream infile;
ofstream outfile;

// Declare big struct of all our data
struct datastruct {string tce,  // TCE string (KIC-PN)
                     comments;  // Comments from ephemeris match
  
                      int kic,  // KIC number
                           pn,  // Planet Number
                  num_planets,  // Number of TCEs for the given KIC
               robo_cent_disp,  // 0 = centroid PC, 1 = centroid FP
             ephem_match_disp;  // 0 = no ephem match, 1 = ephem match identified
               
                double period,  // Period of the system in days from TPS
                        epoch,  // Epoch from TPS
               max_ses_in_mes,  // Max SES actually used in the compuation of the MES
                          mes,  // Max MES from TPS
                     duration,  // Duration of the transit from DV in hours
                       impact,  // Impact parameter of the system from DV
                   dv_sig_pri,  // Significance of primary from  model-shift test on DV flux Data
                   dv_sig_sec,  // Significance of secondary from  model-shift test on DV flux Data
                   dv_sig_ter,  // Significance of tertiary from  model-shift test on DV flux Data
                   dv_sig_pos,  // Significance of positive feature from  model-shift test on DV flux Data
                    dv_sig_fa,  // False Alarm threshold from  model-shift test on DV flux Data
                      dv_fred,  // Red Noise / Gaussian Noise  from  model-shift test on DV flux Data
                    dv_del_fa,  // Delta in sigma value threshold to be considered a distinct feature on DV flux data
                    dv_ph_pri,  // Phase of primary from  model-shift test on DV flux Data
                    dv_ph_sec,  // Phase of secondary from  model-shift test on DV flux Data
                    dv_ph_ter,  // Phase of tertiary from  model-shift test on DV flux Data
                    dv_ph_pos,  // Phase of positive feature from  model-shift test on DV flux Data
              dv_mod_pridepth,  // Depth of primary from DV model-shift
              dv_mod_secdepth,  // Depth of secondary from DV model-shift
                  alt_sig_pri,  // Significance of primary from  model-shift test on Chris' alternate detrended trapezoid fit data
                  alt_sig_sec,  // Significance of secondary from  model-shift test on Chris' alternate detrended trapezoid fit data
                  alt_sig_ter,  // Significance of tertiary from  model-shift test on Chris' alternate detrended trapezoid fit data
                  alt_sig_pos,  // Significance of positive feature from  model-shift test on Chris' alternate detrended trapezoid fit data
                   alt_sig_fa,  // False Alarm threshold from  model-shift test on Chris' alternate detrended trapezoid fit data
                     alt_fred,  // Red Noise / Gaussian Noise  from  model-shift test on Chris' alternate detrended trapezoid fit data
                   alt_del_fa,  // Delta in sigma value threshold to be considered a distinct feature on alternate detrended trapezoid fit data
                   alt_ph_pri,  // Phase of primary from  model-shift test on Chris' alternate detrended trapezoid fit data
                   alt_ph_sec,  // Phase of secondary from  model-shift test on Chris' alternate detrended trapezoid fit data
                   alt_ph_ter,  // Phase of tertiary from  model-shift test on Chris' alternate detrended trapezoid fit data
                   alt_ph_pos,  // Phase of positive feature from  model-shift test on Chris' alternate detrended trapezoid fit data
             alt_mod_pridepth,  // Depth of primary from  model-shift test on Chris' alternate detrended trapezoid fit data      
             alt_mod_secdepth,  // Depth of secondary from  model-shift test on Chris' alternate detrended trapezoid fit data
                     dv_oesig,  // Odd-Even test done on the DV data.
                    alt_oesig,  // Odd-Even test done on the alt data.
                      lpp_tps,  // Susan's LPP value based on DV detrending
                     lpp_trap,  // Susan's LPP value based on alt detrending (Chris)
                       dv_alb,  // Albedo from secondary eclipse monte carlo using model-shift test on DV data
                        dv_rp,  // Planet Radius from secondary eclipse monte carlo using model-shift test on DV data
                      alt_alb,  // Albedo from secondary eclipse monte carlo using model-shift test on alternate data
                       alt_rp,  // Planet Radius from secondary eclipse monte carlo using model-shift test on alternate data
                     marshall;  // Marshall metric for calculating if at least three transits are transit-like
                            };
datastruct data[NMAX];
int sig_sec_eclipse[NMAX],not_transit_like[NMAX],planet_occultation[NMAX],centroid_offset[NMAX],period_is_double[NMAX],ephemeris_match[NMAX];  // Final disposition Flags. Made arrays so TCEs in same system can know about each other.
string nexscidisp;  // NEXSCI disposition
double epochthresh; //Epoch matching threshold for sec. determination
int secfound;  // Int to mark if a secondary has been found in the system
double lppsig;  // LPP sigma value, computed based off number of TCEs
string tmpstr1,tmpstr2;  // Temporary strings
double tmpdob1,tmpdob2,tmpdob3,tmpdob4,tmpdob5,tmpdob6;  // Temporary doubles for computations


int main (int argc, char* argv[])  
  {   
  
  // Get Command Line Inputs
  if(argc>1)
    infilename = argv[1];
  else
    {
    cout << "Name of intput file? ";
    cin >> infilename;
    }    
    
  if(argc>2)
    outfilename = argv[2];
  else
    {
    cout << "Name of output file? ";
    cin >> outfilename;
    }

    
  READDATA();  // Read Input Data
  

  // Figure out KIC, PN numbers, and total number of planets in each system, from TCE string
  for(i=0;i<ntces;i++)
    {
    if(data[i].kic != data[i-1].kic)  // If we're on a new system, figure out how many total planets in system for each TCE
      for(j=1;j<=data[i-1].pn;j++)   // Update number of planets for all previous TCEs in the system
        data[i-j].num_planets = data[i-1].pn; 
    data[i].kic = stoi(data[i].tce.substr(1,9));   // KIC
    data[i].pn  = stoi(data[i].tce.substr(11,2));  // Planet Number
    }
  for(j=1;j<=data[i-1].pn;j++)
    data[i-j].num_planets = data[i-1].pn;  // Update number of planets for the last system
    
          
  // Figure out the LPP threshold based on the number of TCEs
  lppsig = sqrt(2)*INVERFC(1.0/20367);  // Fixing to OPS run number of 20,367 for both ops and injected so they use same thresholds

  
  // Okay let the judging begin!
  outfile.open(outfilename.c_str());  // Open Outfile
  outfile << "# 1:TCE  2:NExScI Disposition  3:Not Transit-Like Flag  4:Significant Secondary Flag  5:Centroid Offset Flag  6:Ephemeris Match Flag  7:Minor Descriptive Flags" << endl;
  for(i=0;i<ntces;i++)
    {      
    data[i].comments="";  
    not_transit_like[i]=sig_sec_eclipse[i]=planet_occultation[i]=period_is_double[i]=centroid_offset[i]=ephemeris_match[i]=0;  // Make sure all flags start at 0 - assumed PC until it fails a test
    
    // Let's keep track if we already found a seconary eclipse in the system or not
    if(i==0 || data[i].kic!=data[i-1].kic || not_transit_like[i-1]==0)  // If it's the first TCE we're looking at, or if it's a new system, or if a transit-like TCE was found in the system since we last found a secondary, start looking for a secondary again.
      secfound=0;
    
    // Check to see if TCE is a secondary eclipse
    if(secfound==0)
      ISSEC();
    
    // If not a secondary, check to see if TCE is Transit-Like
    if(not_transit_like[i]==0 && sig_sec_eclipse[i]==0)
      TRANSITLIKE();
    
    // Now if it is Transit-Like, check to see if there is a significant secondary eclipse
    if(not_transit_like[i]==0 && sig_sec_eclipse[i]==0)
      SECECLIPSE();
    
    // Apply robo centroid disposition
    centroid_offset[i] = data[i].robo_cent_disp;
    
    // Apply ephem match disposition
    ephemeris_match[i] = data[i].ephem_match_disp;
      
    // Make final PC/FP determination
    nexscidisp="PC";
    if(not_transit_like[i]==1)
      nexscidisp="FP";
    if(sig_sec_eclipse[i]==1 && planet_occultation[i]==0)  
      nexscidisp="FP";
    if(centroid_offset[i]==1)
      nexscidisp="FP";
    if(ephemeris_match[i]==1)
      nexscidisp="FP";
     
    // Write output
    outfile << data[i].tce << " " << nexscidisp << " " << not_transit_like[i] << " " << sig_sec_eclipse[i] << " " << centroid_offset[i] << " " << ephemeris_match[i] << " " << data[i].comments << endl;
    }
  
  outfile.close();
  }
  

///////////////////////////////////////////////////////////////////////////////////////////////////
  
// Function to read in the data

void READDATA()
  {
  infile.open(infilename.c_str());
  if(infile.fail()==1) // If file doesn't exist, exit with warning
    {
    cout << infilename << " doesn't exist or cannot open file..." << endl;
    exit(0);
    }
    
  infile.ignore(9E9,'\n');  // Ignore header
  
  i=0;
  infile >> data[i].tce;
  while(!infile.eof())
    {
    infile >> data[i].period;
    infile >> data[i].epoch;
    infile >> data[i].duration;
    infile >> data[i].max_ses_in_mes;
    infile >> data[i].mes;
    infile >> data[i].lpp_tps;
    infile >> data[i].lpp_trap;
    infile >> data[i].marshall;
    infile >> data[i].dv_oesig;
    infile >> data[i].alt_oesig;
    infile >> data[i].dv_sig_pri;
    infile >> data[i].dv_sig_sec;
    infile >> data[i].dv_sig_ter;
    infile >> data[i].dv_sig_pos;
    infile >> data[i].dv_fred;
    infile >> data[i].dv_sig_fa;
    infile >> data[i].dv_del_fa;
    infile >> data[i].alt_sig_pri;
    infile >> data[i].alt_sig_sec;
    infile >> data[i].alt_sig_ter;
    infile >> data[i].alt_sig_pos;
    infile >> data[i].alt_fred;
    infile >> data[i].alt_sig_fa;
    infile >> data[i].alt_del_fa;
    infile >> data[i].dv_rp;
    data[i].alt_rp = data[i].dv_rp;
    infile >> data[i].impact;
    infile >> data[i].dv_alb;
    infile >> data[i].dv_mod_pridepth;
    infile >> data[i].dv_mod_secdepth;
    infile >> data[i].dv_ph_sec;
    infile >> data[i].alt_alb;
    infile >> data[i].alt_mod_pridepth;
    infile >> data[i].alt_mod_secdepth;
    infile >> data[i].alt_ph_sec;
    infile >> data[i].robo_cent_disp;
    infile >> data[i].ephem_match_disp;
    i++;
    infile >> data[i].tce;
    }
  infile.close();
  ntces=i;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

// Function to check if TCE is the secondary eclipse of the system

void ISSEC() {

// If a previous TCE in the system had a sec eclipse detected, and this TCE has the same period / diff epoch as it, then this is the sec eclipse
for(j=1;j<data[i].pn;j++)
  {
  COMPPT(data[i].period,data[i-data[i].pn+j].period,data[i].epoch,data[i-data[i].pn+j].epoch);
  tmpdob3 = data[i-data[i].pn+j].period/data[i].period;  // Modify the period ratio so that it's always the other TCE divided by the current one for this test. That way when tmpdob>=2, it always means the current TCE is half the period or less of the previous one.
  
  epochthresh = WIDTHFAC*data[i-data[i].pn+j].duration/24.0;

  if(not_transit_like[i-data[i].pn+j]==0 && (sig_sec_eclipse[i-data[i].pn+j]==1 || period_is_double[i-data[i].pn+j]==1) && tmpdob1 > PSIG_THRESH && ((fabs(tmpdob4) > epochthresh) || (fabs(tmpdob4) < epochthresh && rint(tmpdob3)>=2 )) && ((fabs(tmpdob6) > epochthresh) || (fabs(tmpdob6) < epochthresh && rint(tmpdob3)>=2)) && ((tmpdob4<0 && tmpdob6<0) || (tmpdob4>0 && tmpdob6>0) || rint(tmpdob3)>=2) )  // Either same period and differnt epoch as a previous TCE, or half the period and same epoch. Both end up corresponding to the secondary eclipse.
    {
    not_transit_like[i]=1;
    sig_sec_eclipse[i]=1;
    if(data[i].comments!="")  data[i].comments+="---";
    data[i].comments+="THIS_TCE_IS_A_SEC";
    j=99;  // Only need to trigger once
    secfound=1;  // Note that we found a secondary for this system, so we don't search for more secondaries
    }
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////
 
// Function to check if the TCE is not transit-like
 
void TRANSITLIKE() {
  
// DV LPP Test
// if(data[i].lpp_tps > lppsig*0.000781) // && data[i].period < 50.0)  // Susan OLD value
if(data[i].lpp_tps > 0.00104504238600969 + lppsig*0.000495720001967656)  // NEW value based on fitting gaussian to injections
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="LPP_DV_TOO_HIGH";
  }


// Alt LPP Test
// if(data[i].lpp_trap > lppsig*0.001001) // && data[i].period < 50.0)  // Susan OLD value
if(data[i].lpp_trap > 0.000667164262937681 + lppsig*0.000417055294849554)  // NEW value based on fitting gaussian to injections
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="LPP_ALT_TOO_HIGH";
  }
  
  
// Check Marshall metric
if(data[i].marshall > 10.0)
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="MARSHALL_FAIL";
  }


// Check is primary is significant in DV 
if(data[i].dv_sig_pri/data[i].dv_fred < data[i].dv_sig_fa && data[i].dv_sig_pri > 0)
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="DV_SIG_PRI_OVER_FRED_TOO_LOW";
  }
// Check is primary is significantly greater than tertiary in DV
if(data[i].dv_sig_pri - data[i].dv_sig_ter < data[i].dv_del_fa && data[i].dv_sig_pri > 0 && data[i].dv_sig_ter > 0)  // 0 indicates NULL result
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="DV_SIG_PRI_MINUS_SIG_TER_TOO_LOW";
  }
// Check is primary is significantly greater than positive in DV
if(data[i].dv_sig_pri - data[i].dv_sig_pos < data[i].dv_del_fa && data[i].dv_sig_pri > 0 && data[i].dv_sig_pos > 0)  // 0 indicates NULL result
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="DV_SIG_PRI_MINUS_SIG_POS_TOO_LOW";
  }
  
  
// Check if primary is significant in ALT
if(data[i].alt_sig_pri/data[i].alt_fred < data[i].alt_sig_fa && data[i].alt_sig_pri > 0)
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="ALT_SIG_PRI_OVER_FRED_TOO_LOW";
  }
// Check is primary is significantly greater than tertiary in ALT
if(data[i].alt_sig_pri - data[i].alt_sig_ter < data[i].alt_del_fa && data[i].alt_sig_pri > 0 && data[i].alt_sig_ter > 0)  // 0 indicates NULL result
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="ALT_SIG_PRI_MINUS_SIG_TER_TOO_LOW";
  }
// Check is primary is significantly greater than positive in ALT
if(data[i].alt_sig_pri - data[i].alt_sig_pos < data[i].alt_del_fa && data[i].alt_sig_pri > 0 && data[i].alt_sig_pos > 0)  // 0 indicates NULL result
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="ALT_SIG_PRI_MINUS_SIG_POS_TOO_LOW";
  }


// Check consistency of transits via SES to MES ratio
if(data[i].max_ses_in_mes/data[i].mes > 0.9 && data[i].period > 90)  // Maybe 0.95?  0.99?   1.0 could work. 
  {
  not_transit_like[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="TRANSITS_NOT_CONSISTENT";
  }


// If a previous TCE in the system has the same period as this one, and it was not transit-like, then this one should be too. Also check if previous TCE was transit-like, if this TCE is triggering off the residual of the earlier TCE.
for(j=1;j<data[i].pn;j++)
  {
  COMPPT(data[i].period,data[i-data[i].pn+j].period,data[i].epoch,data[i-data[i].pn+j].epoch);  // Compute diagnostics on the period and epoch matching
  if(tmpdob1 > PSIG_THRESH)  // This TCE matches the period of a previous TCE in the system
    {
    if(not_transit_like[i-data[i].pn+j]==1)  // If the previous TCE was deemed not transit-like, then this TCE should be not transit-like as well
      {
      not_transit_like[i]=1;
      if(data[i].comments!="")  data[i].comments+="---";
      data[i].comments+="SAME_P_AS_PREV_NTL_TCE";
      j=99;  // Only need to trigger once
      }
    else
      if(fabs(tmpdob4) < WIDTHFAC*data[i-data[i].pn+j].duration/24.0 || fabs(tmpdob6) < WIDTHFAC*data[i-data[i].pn+j].duration/24.0 || (tmpdob4<0 && tmpdob6>0) || (tmpdob4>0 && tmpdob6<0))  // If previous TCE was transit-like, check to see if this TCE is triggering on its residuals, i.e., it's epoch is within two transit durations.
        {
        not_transit_like[i]=1;
        if(data[i].comments!="")  data[i].comments+="---";
        data[i].comments+="RESID_OF_PREV_TCE";
        j=99;  // Only need to trigger once
        }
    }
  }
  
}


///////////////////////////////////////////////////////////////////////////////////////////////////

// Function to check if the TCE has a visible secondary eclipse
  
void SECECLIPSE() {

// Look for secondary in DV detrending  
if(data[i].dv_sig_sec/data[i].dv_fred > data[i].dv_sig_fa && data[i].dv_sig_sec > 0)  // See if secondary is significant
  if(data[i].dv_sig_sec - data[i].dv_sig_ter > data[i].dv_del_fa || data[i].dv_sig_ter <= 0)  // If ter measurement exists, check if sec is more significant
    if(data[i].dv_sig_sec - data[i].dv_sig_pos > data[i].dv_del_fa || data[i].dv_sig_pos <= 0)  // If pos measurement exists, check if sec is more significant
      {
      sig_sec_eclipse[i]=1;
      if(data[i].comments!="")  data[i].comments+="---";
      data[i].comments+="SIG_SEC_IN_DV_MODEL_SHIFT";
      if(data[i].dv_alb > 0.0 && data[i].dv_alb < 1.0 && data[i].dv_rp > 0.0 && data[i].dv_rp < 30.0 && data[i].dv_mod_secdepth < 0.10*data[i].dv_mod_pridepth && data[i].impact < 0.95)  // Check to see if occultation could be due to planet. Only apply to things with less than 30 earth radii, and if the secodary is less than 10% the depth of the primary, and if impact parameter is less than 0.9.
        {
        planet_occultation[i]=1;
        if(data[i].comments!="")  data[i].comments+="---";
        data[i].comments+="DV_SEC_COULD_BE_DUE_TO_PLANET";
        }
      if(fabs(0.5 - data[i].dv_ph_sec)*data[i].period < 0.25*data[i].duration/24.0 && fabs(data[i].dv_sig_pri - data[i].dv_sig_sec) < data[i].dv_del_fa)  // Check to see if secondary could be identical to the transit so that it's really a PC phased at twice the period
        {
        period_is_double[i]=1;
        if(data[i].comments!="")  data[i].comments+="---";
        data[i].comments+="DV_SEC_SAME_DEPTH_AS_PRI_COULD_BE_TWICE_TRUE_PERIOD";
        }
      }


// Look for secondary in Alt detrending
if(data[i].alt_sig_sec/data[i].alt_fred > data[i].alt_sig_fa && data[i].alt_sig_sec > 0)  // See if secondary is significant
  if(data[i].alt_sig_sec - data[i].alt_sig_ter > data[i].alt_del_fa || data[i].alt_sig_ter <= 0)  // If ter measurement exists, check if sec is more significant
    if(data[i].alt_sig_sec - data[i].alt_sig_pos > data[i].alt_del_fa || data[i].alt_sig_pos <= 0)  // If pos measurement exists, check if sec is more significant
      {
      sig_sec_eclipse[i]=1;
      if(data[i].comments!="")  data[i].comments+="---";
      data[i].comments+="SIG_SEC_IN_ALT_MODEL_SHIFT";
      if(data[i].alt_alb > 0.0 && data[i].alt_alb < 1.0 && data[i].alt_rp > 0.0 && data[i].alt_rp < 30.0 && data[i].alt_mod_secdepth < 0.1*data[i].alt_mod_pridepth && data[i].impact < 0.95)  // Check to see if occultation could be due to planet
        {
        planet_occultation[i]=1;
        if(data[i].comments!="")  data[i].comments+="---";
        data[i].comments+="ALT_SEC_COULD_BE_DUE_TO_PLANET";
        }
      if(fabs(0.5 - data[i].alt_ph_sec)*data[i].period < 0.25*data[i].duration/24.0 && fabs(data[i].alt_sig_pri - data[i].alt_sig_sec) < data[i].alt_del_fa)  // Check to see if secondary could be identical to the transit so that it's really a PC phased at twice the period
        {
        period_is_double[i]=1;
        if(data[i].comments!="")  data[i].comments+="---";
        data[i].comments+="ALT_SEC_SAME_DEPTH_AS_PRI_COULD_BE_TWICE_TRUE_PERIOD";
        }
      }


//  Odd-Even Test from DV detrending
if(data[i].period < 90.0 && data[i].dv_oesig > 1.70)
  {
  sig_sec_eclipse[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="DV_ROBO_ODD_EVEN_TEST_FAIL";
  }

//  Odd-Even Test from Chris detrending
if(data[i].period < 90.0 && data[i].alt_oesig > 1.70)
  {
  sig_sec_eclipse[i]=1;
  if(data[i].comments!="")  data[i].comments+="---";
  data[i].comments+="ALT_ROBO_ODD_EVEN_TEST_FAIL";
  }


// Check if subsequent TCE has same period, indicating a secondary eclipse
if(data[i].kic == data[i+1].kic) // Only run this test if there are subsequent TCEs belonging to the same KIC. Mostly just a precaution for injection systems.
  for(j=1;j<=data[i].num_planets-data[i].pn;j++)  // Start looking at the next TCE in the system, and every subsequent TCE, through the last TCE in the system
    {
    COMPPT(data[i].period,data[i+j].period,data[i].epoch,data[i+j].epoch);
    if(tmpdob1 > PSIG_THRESH && (fabs(tmpdob4) > WIDTHFAC*data[i].duration/24.0 || (fabs(tmpdob4) < WIDTHFAC*data[i].duration/24.0 && rint(tmpdob3)>=2 )) && (fabs(tmpdob6) > WIDTHFAC*data[i].duration/24.0 || (fabs(tmpdob6) < WIDTHFAC*data[i].duration/24.0 && rint(tmpdob3)>=2)) && ((tmpdob4<0 && tmpdob6<0) || (tmpdob4>0 && tmpdob6>0) || rint(tmpdob3)>=2))  // Significant period match, at least 2 transit durations away, and either an insigniicant epoch match or more than a day apart. Put in condition so that if close periods, but some drift, will look for drift across the primary transit
      {
      sig_sec_eclipse[i]=1;
      if(data[i].comments!="")  data[i].comments+="---";
      data[i].comments+="OTHER_TCE_AT_SAME_PERIOD_DIFF_EPOCH";
      j=99;  // Only need to trigger this once
      }
    }


// Now if a secondary was detected, but the period could be double the true period, mark it as not actually having a secondary
if(period_is_double[i]==1)
  sig_sec_eclipse[i]=0;
}
  

///////////////////////////////////////////////////////////////////////////////////////////////////

// Function to compute the match significance of two periods and epochs

void COMPPT(double P1, double P2, double T1, double T2) {
  
if(P1 < P2)
  {
  tmpdob1 = PSIG(P1,P2);  // Period match significance
  tmpdob2 = ESIG(T1,T2,P1);  // Epoch match significance
  tmpdob3 = P2/P1;  // Period Ratio
  tmpdob4 = T1 - T2;  // Difference in epoch
  while(tmpdob4 > 0.5*P1)
    tmpdob4 -= P1;
  while(tmpdob4 < -0.5*P1)
    tmpdob4 += P1;
  tmpdob6 = tmpdob4 + int(MISSIONDUR/P2)*(rint(tmpdob3)*P1-P2);  // Difference in epoch by end of mission

  }
else
  {
  tmpdob1 = PSIG(P2,P1);  // Period match significance
  tmpdob2 = ESIG(T2,T1,P2);  // Epoch match significance
  tmpdob3 = P1/P2;  // Period Ratio
  tmpdob4 = T2 - T1;  // Difference in epoch
  while(tmpdob4 > 0.5*P2)
    tmpdob4 -= P2;
  while(tmpdob4 < -0.5*P2)
    tmpdob4 += P2;
  tmpdob6 = tmpdob4 + int(MISSIONDUR/P1)*(rint(tmpdob3)*P2-P1);  // Difference in epoch by end of mission
  }

}

double PSIG(double P1, double P2) {
  return sqrt(2)*INVERFC(fabs((P1-P2)/P1 - rint((P1-P2)/P1)));
}

double ESIG(double E1, double E2, double P) {
  return sqrt(2)*INVERFC(fabs((E1-E2)/P - rint((E1-E2)/P)));
}


///////////////////////////////////////////////////////////////////////////////////////////////////

// The error function 

double INVERFC(double p) {
  double x, err, t, pp;

  if (p >= 2.) return -100.;
  if (p <= 0.0) return 100.;

  pp=(p < 1.0)? p:2.-p;
  t=sqrt(-2.*log(pp/2));

  x= -0.70711*((2.30753+t*0.27061)/(1+t*(0.99229+t*0.04481))-t);

  for (int j=0;j<2;j++) {
    err=erfc(x)-pp;
    x += err/(1.12837916709551257*exp(-x*x)-x*err);
  }
  return (p<1.0 ? x:-x);
}
