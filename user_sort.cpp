/******************************************************************************/
/*** SORTING CODE FOR 184W(alpha,p)187Re                                    ***/
/*** LAST UPDATE: DEC 8, 2015, CECILIE                                      ***/
/******************************************************************************/

using namespace std;

#include "Event.h"
#include "Histogram1D.h"
#include "Histogram2D.h"
#include "IOPrintf.h"
#include "OfflineSorting.h"
#include "Parameters.h"
#include "ParticleRange.h"
#include "SiriusRoutine.h"
#include <TH1D.h>
#include <TH1.h>
#include <TFile.h>


#include "TCutG.h"
#include "TClass.h"
#include "TSystem.h"
#include <TH2.h>
#include <TTree.h>


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <cstring>


#define NDEBUG 1
#include "debug.h"



//! must be 0 to disable making DE:E matrices for each strip of all detectors
// Makes plots of type m_e_de_b0f0...b7f7
#define MAKE_INDIVIDUAL_E_DE_PLOTS 1

//! must be 0 to disable making time evolution plots
// This is to check, for example, if the elastic peak or a gamma peak 
// (which is more likely) moves during the experiment
#define MAKE_TIME_EVOLUTION_PLOTS 1

//! must be 0 to disable making CACTUS time-energy plots
// Makes plots of type m_nai_e_t_0...27
#define MAKE_CACTUS_TIME_ENERGY_PLOTS 1

//! must be 0 to disable checking particles
// If this is on, it checks that 
#define APPLY_PARTICLE_GATE 1




// ########################################################################

//! User sorting routine.
class UserXY : public SiriusRoutine {
public:
    UserXY();
    bool Sort(const Event& event);
    void CreateSpectra();
    bool Command(const std::string& cmd);
    TCutG* get_cut(char* cuts_file_name, char* cut_variable_name);


    
private:
    Histogram2Dp m_back, m_front, m_e_de_strip[8], m_e_de, m_e_de_thick;
#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0)
    Histogram2Dp m_e_de_individual[8][8];
    Histogram1Dp h_ede_individual[8][8];
    // JEM addition 20160229: Make copies of the 2D histograms subject to TCut conditions 
    // Histogram2Dp m_e_de_individual_cut[8][8]; // COMMENTED OUT 20160301: Cannot get it to run with this, instead just applying the cut directly on the m_e_de matrices.
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */
    Histogram2Dp m_nai_t, m_nai_e;
    Histogram2Dp m_nai_e_gate1; // Added by JEM 20160316 to gate on excitation peaks
    Histogram2Dp m_nai_e_gate2, m_nai_e_gate3, m_nai_e_gate4, m_nai_e_gate5; // Added by JEM 20160329 to have more different gates
    Histogram2Dp m_alfna, m_alfna_f0, m_alfna_bg, m_alflabr, m_alflabr_bg;
    
    
    Histogram1Dp h_na_n, h_thick, h_ede, h_ede_r[8], h_ex_r[8], h_de_n, h_e_n;
    Histogram1Dp h_ex;

    Histogram1Dp h_nai_e; // Added by JEM 20160310 to combine all NaI channels
    Histogram1Dp h_nai_e_individual[32]; // Added by JEM 20160311 to plot individual NaI spectra
    

#if defined(MAKE_CACTUS_TIME_ENERGY_PLOTS) && (MAKE_CACTUS_TIME_ENERGY_PLOTS>0)
    Histogram2Dp m_nai_e_t[32], m_nai_e_t_all, m_nai_e_t_c, m_siri_e_t[8], m_siri_e_t_all, m_nai_e_t_LaBr;
#endif /* MAKE_CACTUS_TIME_ENERGY_PLOTS */

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    Histogram2Dp m_nai_t_evol[32], m_nai_e_evol[32];
    Histogram2Dp m_e_evol[8], m_de_evol[8][8], m_ede_evol[8], m_ex_evol;
#endif /* MAKE_TIME_EVOLUTION_PLOTS */

private:
    //! Correction of CACTUS time for CACTUS energy.
    Parameter tnai_corr_enai;
    //! Correction of CACTUS time for SiRi back detector energy.
    Parameter tnai_corr_esi;
    //! Polynomials to calculate excitation energy from SiRi energy (E+DE).
    /*! Contains 8*3 coefficients. */
    Parameter ex_from_ede;
    //! Polynomials to make an empirical correction to the calculated excitation energy.
    /*! Contains 8*2 coefficients. */
    Parameter ex_corr_exp;
    //! Two rectangles to cut away SiRi noise/electrons.
    /*! Contains E-minimum 1, DE-minimum 1, E-minimum 2, DE-minimum 2. */
    Parameter ede_rect;
    //! Thickness centroid and second-order polynomial for SiRi E dependent thickness gate width.
    /*! Contains centroid, constant width, E-scaled width. */
    Parameter thick_range;
    //! The particle range data from zrange.
    ParticleRange particlerange;

    //! Apply energy corrections to CACTUS time.
    /*! \return corrected CACTUS time. */
    double tNaI(double t,    /*!< Uncorrected CACTUS time. */
               double Enai, /*!< Calibrated CACTUS energy in keV. */
               double Esi   /*!< Calibrated SiRi back energy in keV. */);

    //! Apply energy corrections to LaBr time, but only the start signal from the E counter in SiRi.
    /*! There are already CFD corrections on the LaBr3 signals.
        \return corrected CACTUS time. */
    double tLaBr(double t,    /*!< Uncorrected LaBr3 time. */
               double Esi   /*!< Calibrated SiRi back energy in keV. */);
    
    double range(double E /*!< particle energy in keV */)
        { return particlerange.GetRange( (int)E ); }




    // The following added by JEM 20160229 to apply TCuts to matrices:
    // // cuts_file_name, cut_variable_name_p, cut_variable_name_t;

    // TFile *cuts_file = new TFile(cuts_file_name,"read");
    // TFile *cuts_file = new TFile("ipynb/cuts_p_and_t.root","read");


    // TCutG* get_cut(const char* cuts_file_name, const char* cut_variable_name);

    // char* cuts_file_name = "ipynb_test/cuts_p_and_t.root";
    // char* cut_variable_name_p = "cutp2";
    // char* cut_variable_name_t = "cutt1";
    bool apply_TCut = false;
    // TCutG cutp = *get_cut("ipynb_test/cuts_p_and_t.root", "cutp2");
    // TCutG cutt = *get_cut("ipynb_test/cuts_p_and_t.root", "cutt1");
    TCutG cutp = *get_cut("ipynb_test/cutp_20160311.root", "cutp"); // New versions of cuts 
    TCutG cutt = *get_cut("ipynb_test/cutt_20160311.root", "cutt"); // added 20160311.
};

// ########################################################################

bool UserXY::Command(const std::string& cmd)
{
    std::istringstream icmd(cmd.c_str());

    std::string name;
    icmd >> name;
    
    if( name == "rangefile" ) {
        std::string filename;
        icmd >> filename;
        particlerange.Read( filename );
        return true;
    }
    return SiriusRoutine::Command(cmd);
}


// Following added by JEM 20160229 to apply TCut to matrices
// TCut objects can be made and saved in an interactive ROOT session:
// After drawing some 2D plot, click View and then Toolbar. In the 
// toolbar click on the scissors symbol, and click around in the plot where you 
// want to define your cut. Double click to close the cut.
// Then you can save the cut by right-clicking on it and choosing Save As,
// giving it a name that ends with .root to make it a ROOT file.
TCutG* UserXY::get_cut(char* cuts_file_name, char* cut_variable_name)
{
    TFile *cuts_file = new TFile(cuts_file_name,"read");
    TCutG* cut = (TCutG*)cuts_file->Get(cut_variable_name)->Clone();
    // TCutG* cutt = (TCutG*)cuts_file->Get("cutt1")->Clone();
    cuts_file->Close();
    return cut;
}

// ########################################################################

UserXY::UserXY()
    : tnai_corr_enai ( GetParameters(), "tnai_corr_enai", 4 )
    , tnai_corr_esi  ( GetParameters(), "tnai_corr_esi", 4 )
    , ex_from_ede    ( GetParameters(), "ex_from_ede", 8*3 )
    , ex_corr_exp    ( GetParameters(), "ex_corr_exp", 8*2 )
    , ede_rect       ( GetParameters(), "ede_rect", 4 )
    , thick_range    ( GetParameters(), "thick_range", 3 )
{
    ede_rect.Set( "500 250 30 504" ); //"500 250 30 500"
    thick_range.Set( "130  13 0" );
}


// ########################################################################

void UserXY::CreateSpectra()
{
    const int max_e = 40000, max_de = 14000;
    m_back = Mat( "m_back", "back detector energies",
                  2000, 0, max_e, "E(Si) [keV]", 8, 0, 8, "detector nr." );
    m_front = Mat( "m_front", "front detector energies",
                   2000, 0, max_de, "#DeltaE(Si) [keV]", 64, 0, 64, "detector nr." );

#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0)
    // try to change bin counts if not enough statistics for calibration
    for(int b=0; b<8; ++b ) {
        for(int f=0; f<8; ++f ) {
            m_e_de_individual[b][f] = 
                Mat( ioprintf("m_e_de_b%df%d", b, f), ioprintf("#DeltaE : E detector %d strip %d", b, f),
                     2000, 0, max_e, "E(Si) [keV]", 2000, 0, max_de, "#DeltaE(Si) [keV]" );
            h_ede_individual[b][f] =
                Spec( ioprintf("h_ede_b%df%d", b, f), ioprintf("E+#DeltaE detector %d strip %d", b, f),
                      2000, 0, max_e, "E+#DeltaE [keV]" );
            // Added by JEM 20160229 to make 2D histogram versions subject to cuts:
            // m_e_de_individual_cut[b][f] =
            //     Mat( ioprintf("m_e_de_b%df%d", b, f), ioprintf("#DeltaE : E detector %d strip %d, subject to TCut", b, f),
            //          2000, 0, max_e, "E(Si) [keV]", 2000, 0, max_de, "#DeltaE(Si) [keV]" );
        }
    }
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */
    for(int f=0; f<8; ++f ) {
        m_e_de_strip[f] = Mat( ioprintf("m_e_de_f%d", f), ioprintf("E(NaI) : E(Si) strip %d", f),
                               2000, 0, max_e, "E(Si) [keV]", 2000, 0, max_de, "#DeltaE(Si) [keV]" );
    }

    m_e_de = Mat( "m_e_de", "#DeltaE : E for all detectors together",
                  500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    m_e_de_thick = Mat( "m_e_de_thick", "#DeltaE : E for all detectors together, gated on thickness",
                        1000, 0, max_e, "E(Si) [keV]", 1000, 0, max_de, "#DeltaE(Si) [keV]" );
    
    m_nai_t = Mat( "m_nai_t", "t(NaI) matrix", 500, 0,  500, "? [a.u.]",     32,0,32, "det. id.");
    m_nai_e = Mat( "m_nai_e", "E(NaI) matrix", 2000, 0, 14000, "E(NaI) [keV]", 32,0,32, "det. id.");
    m_nai_e_gate1 = Mat( "m_nai_e_gate1", "E(NaI) matrix, Ex gated 7800-8700 keV", 2000, 0, 14000, "E(NaI) [keV]", 32,0,32, "det. id."); // Added by JEM 20160316 to gate on peaks
    m_nai_e_gate2 = Mat( "m_nai_e_gate2", "E(NaI) matrix, Ex gated 8800-10 000 keV", 2000, 0, 14000, "E(NaI) [keV]", 32,0,32, "det. id."); // Added by JEM 20160329
    m_nai_e_gate3 = Mat( "m_nai_e_gate3", "E(NaI) matrix, Ex gated NA-NA keV", 2000, 0, 14000, "E(NaI) [keV]", 32,0,32, "det. id."); // Added by JEM 20160329
    m_nai_e_gate4 = Mat( "m_nai_e_gate4", "E(NaI) matrix, Ex gated NA-NA keV", 2000, 0, 14000, "E(NaI) [keV]", 32,0,32, "det. id."); // Added by JEM 20160329

    h_nai_e = Spec("h_nai_e", "NaI energy", 2000, 0, 14000, "E(NaI) [keV]");// Added by JEM 20160310 to sum all NaI energy spectra
    for (int n=0; n<32; ++n)
    {
        h_nai_e_individual[n] = Spec(ioprintf("h_nai_e_%d", n), ioprintf("NaI energy detector %d", n),
                          2000, 0, 14000, "E(NaI) [keV]");  // Added by JEM 20160311
    }

    m_alfna = Mat( "m_alfna", "E(NaI) : E_{x}",
                   1600, 0, 16000, "E(NaI) [keV]", 500, 0, 20000, "E_{x} [keV]" );

    
    m_alfna_bg = Mat( "m_alfna_bg", "E(NaI) : E_{x} background",
                   1600, 0, 16000, "E(NaI) [keV]", 500, 0, 20000, "E_{x} [keV]" );


    h_na_n = Spec("h_na_n", "NaI multiplicity", 32, 0, 32, "multiplicity");

    h_thick = Spec("h_thick", "apparent #DeltaE thickness", 500, 0, 500, "#DeltaE 'thickness' [um]");
    h_de_n = Spec("h_de_n", "#DeltaE multiplicity", 64, 0, 64, "multiplicity");
    h_e_n = Spec("h_e_n", "E multiplicity", 10, 0, 10, "multiplicity");
    

    for(int f=0; f<8; ++f ) {
        h_ede_r[f] = Spec(ioprintf("h_ede_f%d", f), ioprintf("E+#DelatZE ring %d", f),
                          1000, 0, max_e, "E+#DeltaE [keV]");
        //h_ede_r[f]->SetLineColor(f+1);

        h_ex_r[f] = Spec(ioprintf("h_ex_f%d", f), ioprintf("E_{x} ring %d", f),
                          500, -500, 15000, "E_{x} [keV]");
        //h_ex_r[f]->SetLineColor(f+1);
    }
    h_ede = Spec("h_ede", "E+#DeltaE all detectors", 1000, 0, max_e, "E+#DeltaE [keV]");
    h_ex  = Spec("h_ex", "E_{x} all detectors", 500, -500, 15000, "E_{x} [keV]");

#if defined(MAKE_CACTUS_TIME_ENERGY_PLOTS) && (MAKE_CACTUS_TIME_ENERGY_PLOTS>0)
    const int max_enai = 14000;
    for(int n=0; n<32; ++n ) {
        m_nai_e_t[n] = Mat( ioprintf("m_nai_e_t_%02d", n), ioprintf("t : E NaI %d", n),
                            2000, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );
    }
    m_nai_e_t_all = Mat( "m_nai_e_t", "t : E NaI all together",
                         2000, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );
    m_nai_e_t_LaBr = Mat( "m_nai_e_t_labr", "t : E LaBr all together",
                        2000, 0, max_enai, "E(LaBr) [keV]", 500, 0, 500, "t(LaBr) [a.u.]" );
    m_nai_e_t_c   = Mat( "m_nai_e_t_c", "t : E NaI all together, corrected",
                         2000, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );

    for(int n=0; n<8; ++n ) {
        m_siri_e_t[n] = Mat( ioprintf("m_siri_e_t_b%d", n), ioprintf("t(NaI) : E(Si) detector %d", n),
                             500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(NaI) corr. [a.u.]" );
    }
    m_siri_e_t_all = Mat( "m_siri_e_t", "t(NaI) : E(Si) all detectors",
                          500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(NaI) corr. [a.u.]" );
#endif /* MAKE_CACTUS_TIME_ENERGY_PLOTS */

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    // time evolution plots
    const int MT = 4*24*3600;
    const int NT = 4*24;
    for(int n=0; n<32; ++n ) {
        m_nai_t_evol[n] = Mat( ioprintf("m_nai_t_evol_%02d", n), ioprintf("time : t NaI %d", n),
                            500, 0, 500, "t(NaI) [a.u.]", NT, 0, MT, "wall clock time [s]" );
        m_nai_e_evol[n] = Mat( ioprintf("m_nai_e_evol_%02d", n), ioprintf("time : e NaI %d", n),
                            500, -1000, max_enai-1000, "e(NaI) [keV]", NT, 0, MT, "wall clock time [s]" );
    }
    
    for(int b=0; b<8; ++b ) {
        m_e_evol[b] = Mat( ioprintf("m_e_evol_b%d", b), ioprintf("time : E detector %d", b),
                           500, 0, max_e, "E(Si) [keV]", NT, 0, MT, "wall clock time [s]" );
        for(int f=0; f<8; ++f ) {
            m_de_evol[b][f] = 
                Mat( ioprintf("m_de_evol_b%df%d", b, f), ioprintf("time : #DeltaE detector %d strip %d", b, f),
                     500, 0, max_de, "#DeltaE(Si) [keV]", NT, 0, MT, "wall clock time [s]" );
        }
        m_ede_evol[b] = Mat( ioprintf("m_ede_evol_f%d", b), ioprintf("time : E+#DeltaE ring %d", b),
                             500, 0, max_e, "E+#DeltaE(Si) [keV]", NT, 0, MT, "wall clock time [s]" );
    }
    m_ex_evol  = Mat("m_ex_evol", "time : E_{x} all detectors", 800, -2000, 14000, "E_{x} [keV]",
                     NT, 0, MT, "wall clock time [s]" );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */
}

// ########################################################################

static double _rando = 0;
static double calib(unsigned int raw, double gain, double shift)
{
    return shift + (raw+_rando) * gain;
}

// ########################################################################

double UserXY::tNaI(double t, double Enai, double Esi)
{
    const double c = tnai_corr_enai[0] + tnai_corr_enai[1]/(Enai+tnai_corr_enai[2]) + tnai_corr_enai[3]*Enai;
    const double d = tnai_corr_esi [0] + tnai_corr_esi [1]/(Esi +tnai_corr_esi [2]) + tnai_corr_esi [3]*Esi;
    return t - c - d;
}

// ########################################################################


bool UserXY::Sort(const Event& event)
{
    _rando = drand48() - 0.5;

    
    // ......................................................//
    // MY MODIFIED SORTING ROUTINE
    // GIVES MUCH LESS "PILE-UP" TO THE LEFT OF ELASTIC PEAK
    // ......................................................//
    
    
    const int e_mini = 50;    // E threshold
    int si_e_raw[8];   // E energies
    int raw;
    
    int id_b;
    // ..................................................
    // E DETECTORS
    // Check E counters first, these give the master gate. 
    // Allow for only one E detector within the time range
    // (although there are many events with two and more!)
    // ..................................................
    
    for( int i=0; i<event.n_e; i++ ) {
        //if(event.n_e>1) // only one E with signal, else jump to next event
          //  return true;
        
        int id = event.e[i].chn;
        if( !(id&1) || id>= 16 )
            continue; // ignore guard rings. Guard rings in id 0,2,4,6,8,10,12,14, detectors in id 1,3,5,7,9,11,13,15
        
        id >>= 1; // detector id is now transformed into detector number 0..7
        id_b = id;  // keep the E ID
        
        // only keep raw E here, we don't know which front detector
        // has fired, so we also don't know which coefficients to use
        // to calibrate the back detector
        raw = event.e[i].adc;
        if( raw >= e_mini ) // check against E threshold
            si_e_raw[id] = raw; // if above, OK
        else{
            si_e_raw[id] = 0;   // if below, set energy to 0 and continue to next E counter
            continue;
        }
        
    }
    // approximate calibration
    m_back->Fill( (int)calib( raw, gain_e[8*id_b], shift_e[8*id_b] ), id_b );
    h_e_n->Fill(event.n_e);
    
    // ..................................................
    // DELTA E DETECTORS. 
    // NOW WE KNOW WHICH BACK DETECTOR HAS FIRED, SO WE
    // CAN CHECK THAT WE HAVE THE CORRECT FRONT DETECTOR.
    // CORRELATIONS:
    // id_back = 0, id_front = 0...7
    // id_back = 1, id_front = 8...15
    // id_back = 2, id_front = 16...23
    // id_back = 3, id_front = 24...31
    // id_back = 4, id_front = 32...39
    // id_back = 5, id_front = 40...47
    // id_back = 6, id_front = 48...55
    // id_back = 7, id_front = 56...63
    // ..................................................
    
    int si_goodcount = 0, dei=-1, ei=-1;
    double de = 0;
    int i_front_start = id_b*8;
    int i_front_stop  = (id_b*8) + 7;
    
    //    std::cout << " ............................................." << std::endl;
    //    std::cout << " Back detector fired: " << id_b <<  std::endl;
    //    std::cout << " Start front: " << i_front_start << ", stop front: " << i_front_stop << std::endl;
    //    std::cout << " Number of Delta E events:" << event.n_de << std::endl;
    for(int i=0;i<event.n_de;i++){
        const int id   = event.de[i].chn;
        //        std::cout << " Front ID fired: " << id << std::endl;
        
        if(id<i_front_start || id>i_front_stop) // check id of front against the id of back - must match
            continue;
        
        //        std::cout << " Back ID: " << id_b << ", front strip " << id << ", energy front:" << raw << std::endl;
        
        const unsigned int raw = event.de[i].adc;
        const double de_cal = calib( raw, gain_de[id], shift_de[id] );
        if(de_cal < 50) // threshold on delta E's
            continue;
        
        m_front->Fill( (int)de_cal, id );   
        
        const int id_f = id % 8; // from here, the front strip has ID 0..7
        
        // Now ensure that only one strip has fired
        if(si_goodcount < 2){
            ei  = id_b;
            dei = id_f;
            de  = de_cal;
            si_goodcount += 1;
        }
    }
    
    
    h_de_n->Fill(event.n_de);  

    
    if( APPLY_PARTICLE_GATE && si_goodcount != 1 )
        // no detector above threshold, reject event
        return true;
 

    const double e  = calib( si_e_raw[ei], gain_e[8*ei+dei], shift_e[8*ei+dei] );
    const int e_int = int(e), de_int = int(de);

    // ..................................................

    // make DE:E matrices
#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0) 
    // m_e_de_individual[ei][dei]->Fill( e_int, de_int );
    // The following is added by JEM 20160229: Filling matrices subject to TCut:
    if (!apply_TCut)
    {
        m_e_de_individual[ei][dei]->Fill( e_int, de_int );
    }
    else if (apply_TCut && (cutp.IsInside(e_int, de_int) || cutt.IsInside(e_int, de_int)))
    {
        // cout << "YES" << endl;
        m_e_de_individual[ei][dei]->Fill( e_int, de_int );
        // m_e_de_individual_cut[ei][dei]->Fill( e_int, de_int ); // COMMENTED OUT 20160301: Cannot get it to run with two sets of matrices, instead applying the cut directly to the m_e_de_individual and _strip
    }
    // else
    // {
    //     cout << "NO" << endl;
    // }
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */
    if (!apply_TCut)
    {
        m_e_de_strip[dei]->Fill( e_int, de_int );
        m_e_de->Fill( e_int, de_int );
    }
    else if (apply_TCut && (cutp.IsInside(e_int, de_int) || cutt.IsInside(e_int, de_int))) // If test added by JEM 20160301 to apply cut.
    {
        m_e_de_strip[dei]->Fill( e_int, de_int );
        m_e_de->Fill( e_int, de_int );
    }
    //const double thick = range(e+de)-range(e);
    const double thick = range(e+de)-range(e);
    h_thick->Fill( (int)thick );

    const double thick_dev = thick_range[1] + thick_range[2]*e;
    const bool have_pp = fabs(thick-thick_range[0])<thick_dev;
    if( APPLY_PARTICLE_GATE && !have_pp )
        return true;

    m_e_de_thick->Fill( e_int, de_int );
    double ede = e+de;
    int   ede_int = (int)ede;
    
    // if f3 is at 42 degrees instead of 46 deg, target shifted 5 mm backwards with respect to SiRi
    if(dei==3){
        ede = ede+75;
        ede_int = (int)ede;
    }

    h_ede->Fill( ede_int );
    h_ede_r[dei]->Fill( ede_int );
#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0)
    h_ede_individual[ei][dei]->Fill( ede_int );
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */

    // fit of kinz Ex(E+DE)
    const double ex_theo = ex_from_ede[3*dei+0] + (ede)*(ex_from_ede[3*dei+1] + (ede)*ex_from_ede[3*dei+2]);
    //const double ex_theo = ex_from_ede.Poly(ede, 3*dei, 3);

    // make experimental corrections
    const double ex = ex_corr_exp[2*dei]+ex_corr_exp[2*dei+1]*ex_theo;
    const int   ex_int = (int)ex;

    h_ex->Fill( ex_int );
    h_ex_r[dei]->Fill( ex_int );

    // ..................................................

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    const int timediff = Timediff(event);
    m_ex_evol->Fill( ex_int, timediff );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */

    // ..................................................
    
    
    /*** Now start with the gamma detectors ***/
    h_na_n->Fill(event.n_na);
    
    // Vectors to keep the NaI ID, energy, and time
    int nai_id[32] = {-1};
    int nai_e[32] = {-1};
    int nai_t[32] = {-1};
    int numbers_of_nais_with_signal = 0;

    for( int i=0; i<event.n_na; i++ ) {
        const int id = event.na[i].chn;
        if( event.na[i].adc <= 0 )
            continue;
        //std::cout << id << std::endl;
        double na_e = calib( (int)event.na[i].adc, gain_na[id], shift_na[id] );
        //cout << " Detector ID = " << id << ", energy = " << na_e << " keV." << endl;
        
        
        const int   na_e_int = (int)na_e;
        




        if( event.na[i].tdc <= 0 )
            continue;
        const double na_t = calib( (int)event.na[i].tdc/8, gain_tna[id], shift_tna[id] ); 
        const int   na_t_int = (int)na_t;
        const int   na_t_c = (int)tNaI(na_t, na_e, e);
        if(id >= 0 && id<31 ) // && ex_int>8705 && ex_int<9002
            m_nai_t->Fill( na_t_c, id );

        // Keep id, energy, and time of those NaI's that gave an energy and time signal
        for(int j=0;j<31;j++){
            if(id==j && na_e>0 && na_t_c>0){  // ADC ID 0-31
                nai_id[j] = id;
                nai_e[j] = (int) na_e;
                nai_t[j] = na_t_c;
            }
        }


        // === NaI energy spectra: ===
        // Everything:
        m_nai_e->Fill( na_e_int, id );
        // Gate on peak at ~8 MeV which we believe to be contamination, maybe the nice 5 MeV line?
        if (ex_int > 7800 && ex_int < 8700 && na_t_c>135 && na_t_c<230) // 20160331: Testing time gates
            m_nai_e_gate1->Fill( na_e_int, id );
        else if(ex_int > 8800 && ex_int < 10000 && na_t_c>135 && na_t_c<230) // Added 20160329. Have also allocated two more gate matrices which may be filled just below.
            m_nai_e_gate2->Fill( na_e_int, id );
        else if(ex_int > 6700 && ex_int < 7800 && na_t_c>135 && na_t_c<230) // Added 20160329. Have also allocated two more gate matrices which may be filled just below.
            m_nai_e_gate3->Fill( na_e_int, id );
        else if(ex_int > 5200 && ex_int < 6700 && na_t_c>135 && na_t_c<230) // Added 20160329. Have also allocated two more gate matrices which may be filled just below.
            m_nai_e_gate4->Fill( na_e_int, id );



        // PUT GATE ON PARTICLE PEAK 
//        if(ex_int>1650 && ex_int<1918)    // 1st Ex in 28Si
        //if(ex_int>8705 && ex_int<9002)      // peak at 8904 keV in 28Si
        //if(ex_int>3550 && ex_int<370



        h_nai_e->Fill( na_e_int ); // Added by JEM 20160310 to combine all NaI energy spectra
        h_nai_e_individual[id]->Fill( na_e_int ); // Added by JEM 20160310 to plot individual NaI spectra


#if defined(MAKE_CACTUS_TIME_ENERGY_PLOTS) && (MAKE_CACTUS_TIME_ENERGY_PLOTS>0)
        //if(ex_int>600 && ex_int<2800)      // gate on particle spectrum
            m_nai_e_t[id] ->Fill( na_e_int,  na_t_int );
//        if(ex_int>8705 && ex_int<9002 && id<26) {     // gate on peak at 8904 keV in 28Si, only NaI
        if(id<31){//ex_int<10000 && 
            m_nai_e_t_all ->Fill( na_e_int,  na_t_int );
            m_nai_e_t_c   ->Fill( na_e_int,  na_t_c );
        }
        m_siri_e_t[ei]->Fill( e_int, na_t_c );
        m_siri_e_t_all->Fill( e_int, na_t_c );
#endif /* MAKE_CACTUS_TIME_ENERGY_PLOTS */

        // ..................................................
       
        /*** HERE COMES THE MAIN MATRIX FOR NaI ***/
        int weight = 1;
        if( id >= 0 && id<31 && na_t_c>135 && na_t_c<230) {// PROMPT TIME GATES 
            m_alfna->Fill( na_e_int, ex_int, 1 );            
        } else if( id >= 0 && id<31 && na_t_c>270 && na_t_c<365) {// RANDOM TIME GATES 
            m_alfna->Fill( na_e_int, ex_int, -1 );
            m_alfna_bg->Fill( na_e_int, ex_int );
            weight = -1;
        }
		

       
 
#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
        m_nai_e_evol[id]->Fill( na_e_int, timediff, 1 );
        m_nai_t_evol[id]->Fill( na_t_c,   timediff );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */
        numbers_of_nais_with_signal += 1;

    } /* END OF NaI LOOP */
    
    
#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    m_e_evol  [ei]     ->Fill( e_int,   timediff );
    m_de_evol [ei][dei]->Fill( de_int,  timediff );
    m_ede_evol[dei]    ->Fill( ede_int, timediff );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */
    


    return true;
}

// ########################################################################
// ########################################################################
// ########################################################################

int main(int argc, char* argv[])
{


    

    return OfflineSorting::Run(new UserXY(), argc, argv );
}



