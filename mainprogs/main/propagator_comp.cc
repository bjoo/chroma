// $Id: propagator_comp.cc,v 1.4 2004-05-11 13:29:29 bjoo Exp $
// $Log: propagator_comp.cc,v $
// Revision 1.4  2004-05-11 13:29:29  bjoo
// Tested lg5eps_w.cc
//
// Revision 1.3  2004/04/28 14:04:17  bjoo
// Fixed propagator_comp for new propagator structure. Propagator stuff appears to work. Yet to check seqsource
//
// Revision 1.2  2004/04/26 11:19:13  bjoo
// Added support for Reading EV-s in SZIN format. Must provide e-values in XML tho.
//
// Revision 1.1  2004/04/23 11:23:38  bjoo
// Added component based propagator (non multishift) and added Zolotarev5D operator to propagator and propagator_comp. Reworked propagator collection scripts
//
// Revision 1.2  2004/04/22 16:49:23  bjoo
// Added overlap_state_info Param struct and gauge startup convenience function. Tidyed up propagator zolo4d case
//
// Revision 1.1  2004/04/20 13:08:12  bjoo
// Added multi mass component based solves and propagator collection program
//
// Revision 1.3  2004/04/16 22:03:59  bjoo
// Added Zolo 4D test files and tidyed up
//
// Revision 1.2  2004/04/16 20:18:03  bjoo
// Zolo seems to work now
//
// Revision 1.1  2004/04/16 17:04:49  bjoo
// Added multi_propagator for Zolo4D multi massery. Seems to work even
//
// Revision 1.49  2004/04/15 14:43:25  bjoo
// Added generalised qprop_io FermAct reading
//
// Revision 1.48  2004/04/06 04:20:33  edwards
// Added SZINQIO support.
//
// Revision 1.47  2004/04/01 18:10:22  edwards
// Added support for non-relativistic quark props.
//
// Revision 1.46  2004/02/23 03:13:58  edwards
// Major overhaul of input/output model! Now using EXCLUSIVELY
// SciDAC propagator format for propagators. Now, Param part of input
// files directly matches source/sink/propagator/seqprop headers
// of propagators. All ``known'' input of a propagator is derived
// from its header(s) and used for subsequent calculations.
//
// Revision 1.45  2004/02/11 12:51:35  bjoo
// Stripped out Read() and Write()
//
// Revision 1.44  2004/02/07 04:51:58  edwards
// Removed SSE hack - pushed it down into the SSE code where it belongs.
//
// Revision 1.43  2004/02/06 22:31:00  edwards
// Put in sse hack for the short term.
//
// Revision 1.42  2004/02/06 17:39:05  edwards
// Added a flush to xml_out.
//
// Revision 1.41  2004/02/05 04:18:56  edwards
// Changed call of quarkProp4 to write to xml_out instead of xml buffer.
//
// Revision 1.40  2004/02/04 17:01:55  edwards
// Changed getSubset() to the now correct getSet().
//
// Revision 1.39  2004/01/31 23:22:01  edwards
// Added proginfo call.
//
// Revision 1.38  2004/01/30 21:35:22  kostas
// added calls to calculate mres for dwf
// 
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>
#include <iomanip>
using namespace std;
#include "chroma.h"

using namespace QDP;


// define MRES_CALCULATION in order to run the code computing the residual mass
// and the pseudoscalar to conserved axial current correlator
#define MRES_CALCULATION

/*
 * Input 
 */
struct Prop_t
{
  string          source_file;
  string          prop_file;
  QDP_volfmt_t    prop_volfmt;
};


struct Component_t { 
  int color;
  int spin;
};

struct PropagatorComponent_input_t
{
  ChromaProp_t     param;
  Cfg_t            cfg;
  Prop_t           prop;
  multi1d<Component_t> components;
};


void read(XMLReader& xml, const string& path, Component_t &comp)
{
  XMLReader top(xml,path);
  try {
    read(top, "color", comp.color);
    read(top, "spin",  comp.spin);
  }
  catch( const string& e ) {
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }  
  if( comp.color < 0 || comp.color >= Nc ) { 
    QDPIO::cerr << "Component color >= Nc. color = " << comp.color << endl;
    QDP_abort(1);
  }

  if( comp.spin < 0 || comp.spin >= Ns ) { 
    QDPIO::cerr << "Component spin >= Ns.  spin = " << comp.spin << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const Component_t &comp)
{
  
  push( xml, path );

  write(xml, "color", comp.color);
  write(xml, "spin",  comp.spin);
  
  pop( xml );
}


// Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "prop_volfmt", input.prop_volfmt);  // singlefile or multifile
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, PropagatorComponent_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the source/propagator info
    read(inputtop, "Prop", input.prop);

    read(inputtop, "Components", input.components);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}

// Forward declaration
void saveComponent(const ChromaProp_t& param, 
		   const Prop_t& prop,
		   XMLReader& source_record_xml,
		   const Component_t& component,
		   XMLReader& gauge_xml,
		   XMLWriter& xml_out,
		   const LatticeFermion& psi,
		   bool make_sourceP,
		   bool seqsourceP);
//! Propagator generation
/*! \defgroup propagator Propagator generation
 *  \ingroup main
 *
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  PropagatorComponent_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  try {
    read(xml_in, "/propagatorComp", input);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << "propagatorComp" << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_file_xml, source_record_xml;

  // ONLY SciDAC mode is supported for propagators!!
  readQprop(source_file_xml, 
	    source_record_xml, quark_prop_source,
	    input.prop.source_file, QDPIO_SERIAL);

  // Try to invert this record XML into a source struct
  // Also pull out the id of this source
  int t0;
  int j_decay;
  multi1d<int> boundary;
  bool make_sourceP = false;;
  bool seqsourceP = false;
  {
    // ONLY SciDAC mode is supported for propagators!!
    QDPIO::cout << "Attempt to read source" << endl;
    readQprop(source_file_xml, 
	      source_record_xml, quark_prop_source,
	      input.prop.source_file, QDPIO_SERIAL);
    QDPIO::cout << "Forward propagator successfully read" << endl;

    // Try to invert this record XML into a source struct
    try
    {
      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	PropSource_t source_header;

	read(source_record_xml, "/MakeSource/PropSource", source_header);
	j_decay = source_header.j_decay;
	t0 = source_header.t_source[j_decay];
	boundary = input.param.boundary;
	make_sourceP = true;
      }
      else if (source_record_xml.count("/SequentialSource") != 0)
      {
	ChromaProp_t prop_header;
	PropSource_t source_header;
	SeqSource_t seqsource_header;

	read(source_record_xml, "/SequentialSource/SeqSource", seqsource_header);
	// Any source header will do for j_decay
	read(source_record_xml, "/SequentialSource/ForwardProps/elem[1]/ForwardProp", 
	     prop_header);
	read(source_record_xml, "/SequentialSource/ForwardProps/elem[1]/PropSource", 
	     source_header);
	j_decay = source_header.j_decay;
	t0 = seqsource_header.t_sink;
	boundary = prop_header.boundary;
	seqsourceP = true;
      }
      else
	throw "No appropriate header found";
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting source_header: " << e << endl;
      QDP_abort(1);
    }
  }    

  // Sanity check
  if (seqsourceP)
  {
    for(int i=0; i < boundary.size(); ++i)
      if (boundary[i] != input.param.boundary[i])
      {
	QDPIO::cerr << "Incompatible boundary between input and seqsource" << endl;
	QDP_abort(1);
      }
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "propagatorComp");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_file_info", source_file_xml);
  write(xml_out, "Source_record_info", source_record_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  // Sanity check - write out the norm2 of the source in the Nd-1 direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);

    multi1d<Double> source_corr = sumMulti(localNorm2(quark_prop_source), 
					   phases.getSet());

    push(xml_out, "Source_correlator");
    write(xml_out, "source_corr", source_corr);
    pop(xml_out);
  }

  xml_out.flush();

  /*
   * Construct fermionic BC. Need one for LatticeFermion and multi1d<LatticeFermion>
   * Note, the handle is on an ABSTRACT type
   */
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(input.param.boundary));
  //
  // Initialize fermion action
  //
  FermionAction<LatticeFermion>* S_f_ptr = 0;
  FermionAction< multi1d<LatticeFermion> >* S_f_a_ptr = 0;

  switch (input.param.FermActHandle->getFermActType() ) {
  case FERM_ACT_WILSON:
    {
      const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_WILSON" << endl;
      S_f_ptr = new EvenOddPrecWilsonFermAct(fbc, wils.Mass,
					     wils.anisoParam);
    }
    break;
    
  case FERM_ACT_UNPRECONDITIONED_WILSON:
    {
      const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_UNPRECONDITIONED_WILSON" << endl;
      S_f_ptr = new UnprecWilsonFermAct(fbc, wils.Mass);
    }
    break;
    
  case FERM_ACT_ZOLOTAREV_4D:
    {
      QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      
      // Construct Fermact -- now uses constructor from the zolo4d params
      // struct
      S_f_ptr = new Zolotarev4DFermAct(fbc, zolo4d, xml_out);
    }
    break;
  
  case FERM_ACT_ZOLOTAREV_5D:
    {
      QDPIO::cout << "FERM_ACT_ZOLOTAREV_5D" << endl;
      const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(input.param.FermActHandle));
    
      // Construct Fermact -- now uses constructor from the zolo4d params
      // struct
      S_f_a_ptr = new Zolotarev5DFermActArray(fbc_a, fbc, zolo5d, xml_out);
    }
    break;

  case FERM_ACT_DWF:
    {
      const DWFFermActParams& dwf = dynamic_cast<const DWFFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_DWF" << endl;
      S_f_a_ptr = new EvenOddPrecDWFermActArray(fbc_a,
						dwf.chiralParam.OverMass, 
						dwf.Mass, 
						dwf.chiralParam.N5);
    }
  break;

  case FERM_ACT_UNPRECONDITIONED_DWF:
    {
      const DWFFermActParams& dwf = dynamic_cast<const DWFFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_UNPRECONDITONED_DWF" << endl;
      S_f_a_ptr = new UnprecDWFermActArray(fbc_a,
					   dwf.chiralParam.OverMass, 
					   dwf.Mass, 
					   dwf.chiralParam.N5);
    }
    break;
  default:
    QDPIO::cerr << "Unsupported fermion action" << endl;
    QDP_abort(1);
  }

  // Create a useable handle on the action
  // The handle now owns the pointer
  Handle< FermionAction<LatticeFermion> > S_f(S_f_ptr);
  Handle< FermionAction< multi1d<LatticeFermion> > > S_f_a(S_f_a_ptr);
  
  // FIrst we have to set up the state -- this is fermact dependent
  const ConnectState *state_ptr;


  switch(input.param.FermActHandle->getFermActType()) {
  case FERM_ACT_WILSON:
    state_ptr = S_f->createState(u);
    break;
  case FERM_ACT_UNPRECONDITIONED_WILSON:
    state_ptr = S_f->createState(u);
    break;
  case FERM_ACT_DWF:
    state_ptr = S_f->createState(u);
    break;
  case FERM_ACT_UNPRECONDITIONED_DWF:
    state_ptr = S_f->createState(u);
    break;
    
  case FERM_ACT_ZOLOTAREV_4D:
    {    
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      const Zolotarev4DFermAct& S_zolo4 = dynamic_cast<const Zolotarev4DFermAct&>(*S_f);
      
      state_ptr = S_zolo4.createState(u, zolo4d.StateInfo, xml_out,zolo4d.AuxFermActHandle->getMass());
      
    }
    break;
  case FERM_ACT_ZOLOTAREV_5D:
    {
      const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(input.param.FermActHandle));
      const Zolotarev5DFermActArray& S_zolo5 = dynamic_cast<const Zolotarev5DFermActArray&>(*S_f_a);
      
      
      state_ptr = S_zolo5.createState(u, zolo5d.StateInfo, xml_out, zolo5d.AuxFermActHandle->getMass());
    }
    break;
    
  default:
    QDPIO::cerr << "Unsupported fermion action (state creation)" << endl;
    QDP_abort(1);
  }
  
  // Now do the quarkprop here again depending on which of the
  // action pointers is not null
  
  Handle<const ConnectState> state(state_ptr);  // inserts any BC


  //
  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //

  
  LatticeFermion psi;
  int ncg_had = 0;
  for(int comp = 0; comp < input.components.size(); comp++) { 
    LatticeFermion chi;

    PropToFerm(quark_prop_source, chi, input.components[comp].color,
	       input.components[comp].spin);
    

    // Zero out initial guess
    psi = zero;
	

    // Normalise source
    Real fact = Real(1) / sqrt(norm2(chi));
    chi *= fact;

    int n_count = 0;
      
    if( S_f_ptr != 0x0 ) { 
      
      S_f->qprop(psi,
		 state,
		 chi,
		 input.param.invParam.invType, 
		 input.param.invParam.RsdCG, 
		 input.param.invParam.MaxCG,
		 n_count);
	
    }
    else if ( S_f_a_ptr != 0x0 ) { 
	
      S_f_a->qprop(psi,
		   state,
		   chi,
		   input.param.invParam.invType, 
		   input.param.invParam.RsdCG, 
		   input.param.invParam.MaxCG,
		   n_count);
    }
    else {
      QDPIO::cerr << "Both S_f_ptr and S_f_a_ptr == 0 " << endl;
      QDP_abort(1);
    }
      
    ncg_had += n_count;
      
    fact = Real(1) / fact;
    psi *= fact;
      
    push(xml_out,"Relaxation_Iterations");
    write(xml_out, "ncg_had", n_count);
    pop(xml_out);
    
      
    saveComponent( input.param,
		   input.prop,
		   source_record_xml,
		   input.components[comp],
		   gauge_xml,
		   xml_out,
		   psi,
		   make_sourceP,
		   seqsourceP);

  }
  
  pop(xml_out);  // propagator
  
  xml_out.close();
  xml_in.close();
  
  // Time to bolt
  QDP_finalize();
  
  exit(0);
}

void saveComponent(const ChromaProp_t& param, 
		   const Prop_t& prop,
		   XMLReader& source_record_xml,
		   const Component_t& component,
		   XMLReader& gauge_xml,
		   XMLWriter& xml_out,
		   const LatticeFermion& psi,
		   bool make_sourceP,
		   bool seqsourceP)
{

  SftMom phases(0, true, Nd-1);
  multi1d<Double> prop_corr = sumMulti(localNorm2(psi), 
				       phases.getSet());
	    
  push(xml_out, "PropComp_correlator");
  write(xml_out, "Mass", param.FermActHandle->getMass());
  write(xml_out, "prop_corr", prop_corr);
  pop(xml_out);

  
    
  XMLBufferWriter file_xml;
  int id = 0;
  push(file_xml, "propagatorComponent");
  // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
  write(file_xml, "id", id);
  pop(file_xml);
  
  XMLBufferWriter record_xml;

  if( make_sourceP) {
    

    XMLReader xml_tmp(source_record_xml, "/MakeSource");

    push(record_xml, "Propagator");    
    write(record_xml, "ForwardProp", param);
    record_xml << xml_tmp;
    write(record_xml, "Component", component);
    pop(record_xml);
  }
  else {
    XMLReader xml_tmp(source_record_xml, "/SequentialSource");

    push(record_xml, "SeqProp");    
    write(record_xml, "ForwardProp", param);
    record_xml << xml_tmp;
    write(record_xml, "Component", component);
    pop(record_xml);
  }
  ostringstream outfile;

  // Zero suffix for compatibility with multi mass version
  outfile << prop.prop_file << "_component_s" << component.spin
	  << "_c" << component.color ;
	  
	  
  QDPIO::cout << "Attempting to write " << outfile.str() << endl;
  
  // Write the source
  writeFermion(file_xml, record_xml, psi,
	       outfile.str(), prop.prop_volfmt, QDPIO_SERIAL);
}
