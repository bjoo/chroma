#include "chromabase.h"

#include "handle.h"

#include "util/gauge/reunit.h"
#include "gtest/gtest.h"

#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/fermacts/eoprec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"



#include "io/xml_group_reader.h"
#include "./qudanef_xml.h"

#include "quda.h"


using namespace Chroma;
using namespace QDP;
using namespace QUDANefTesting;

template<typename TestType>
class FixtureT : public TestType {
public:
  using T = LatticeFermion;
  using Q = multi1d<LatticeColorMatrix>;
  using P = multi1d<LatticeColorMatrix>;

  
  using S_prec_t = EvenOddPrecNEFFermActArray;
  using S_unprec_t = UnprecNEFFermActArray;


  using M_unprec_t = UnprecNEFDWLinOpArray;
  using M_prec_t =  EvenOddPrecNEFDWLinOpArray;
  using Solver_t = LinOpSystemSolverArray<T>;
  
  void SetUp() {
    u.resize(Nd);
    for(int mu=0; mu < Nd; ++mu) {
      gaussian(u[mu]);
      reunit(u[mu]);
    }

    std::istringstream input_fermact(fermact_xml);
    std::istringstream input_fermact_unprec(unprec_fermact_xml);

    QDPIO::cout << "Parsing Prec. Fermact XML" << std::endl;
    XMLReader xml_in_fermact(input_fermact);

    QDPIO::cout << "Parsing Unprec. Fermact XML" << std::endl;
    XMLReader xml_in_fermact_unprec(input_fermact_unprec);


    QDPIO::cout << "Parsing Cuda CG Solver XML" << std::endl;
    std::istringstream input_solver(inv_param_quda_cg_xml);
    xml_in_quda_solver.open(input_solver);

    QDPIO::cout << "Parsing Chroma CG Solver XML" << std::endl;
    std::istringstream input_chroma_solver(inv_param_chroma_cg_xml);
    xml_in_chroma_solver.open(input_chroma_solver);

    QDPIO::cout << "Creating Prec NEF action" << std::endl;
						 
    S_prec = dynamic_cast<S_prec_t*>(TheFermionActionFactory::Instance().createObject("NEF",
										      xml_in_fermact,
										      "FermionAction"));
						 
						 
    QDPIO::cout << "Creating Unprec FermAct" << std::endl;
											    
    S_unprec = dynamic_cast<S_unprec_t*>(TheFermionActionFactory::Instance().createObject("UNPRECONDITIONED_NEF",
											xml_in_fermact_unprec,
											  "FermionAction"));

    state = S_prec->createState(u);

    QDPIO::cout << "Creating Prec LinOp" << std::endl;						
    M_prec = dynamic_cast<M_prec_t *>(S_prec->linOp(state));

    QDPIO::cout << "Creating UnPrec LinOp" << std::endl;						
    M_unprec = dynamic_cast<M_unprec_t*>(S_unprec->linOp(state));
    

  }


  void TearDown() {}

  Q u;

  Handle< FermState<T,P,Q> > state;
  Handle< S_prec_t > S_prec;
  Handle< S_unprec_t > S_unprec;

  Handle< M_prec_t > M_prec;
  Handle <M_unprec_t> M_unprec;

  XMLReader xml_in_quda_solver;
  XMLReader xml_in_chroma_solver;
  
  
};

class TestFixture : public FixtureT<::testing::Test> {};


TEST_F(TestFixture, CheckChromaOpConsistency)
{
  auto N5 = M_prec->size();
  multi1d<T> source(N5);
  multi1d<T> res1(N5);
  multi1d<T> res2(N5);

  for(auto i=0; i < N5; ++i) {
    gaussian(source[i]);
    res1[i] = zero;
    res2[i] = zero;
  }
  QDPIO::cout << "Applying UnprecOp" << std::endl;
  (*M_unprec)(res1, source, PLUS);

  QDPIO::cout << "Applying unprecOp of the PrecOp" << std::endl;
  M_prec->unprecLinOp(res2, source, PLUS);

  Double num_dof = Double(Layout::vol()*4*3*2);
  for(auto i=0; i < N5; ++i) {
    Double diff = sqrt(norm2(res1[i] - res2[i]));
    QDPIO::cout << "i=" << i << " diff = " << diff << " diff/d.o.f = " << diff/num_dof << std::endl;
  } 
}

TEST_F(TestFixture, CheckPrecSolver)
{
  GroupXML_t inv_param = readXMLGroup(xml_in_quda_solver, "//InvertParam", "invType");

  // This instantiates a QUDA solver and sets up the QUDA stuff in the constructor.
  // We can start calling QUDA functions after this.
  
  Handle< Solver_t > Solv_quda = S_prec->invLinOp(state, inv_param);

  auto N5 = M_prec->size();

  multi1d<T> soln(N5);
  multi1d<T> rhs(N5);
  multi1d<T> check(N5);
  for(int i=0; i < N5; ++i) {
    soln[i] = zero;
    check[i] = zero;
    rhs[i] = zero;
    gaussian(rhs[i], rb[1]);
  }

  (*Solv_quda)(soln, rhs);

  (*M_prec)(check, soln, PLUS);

  Double diff_tot(0);
  for(int i=0; i < N5; ++i) {
    Double diff = norm2(check[i]-rhs[i], rb[1]);
    Double diff_b = norm2(rhs[i], rb[1]);
    diff_tot += diff;
    QDPIO::cout << "i=0 Diff = " << sqrt(diff) << " relative= " << sqrt(diff/diff_b) << std::endl;
  }
  QDPIO::cout << "Diff Tot = " << sqrt(diff_tot) << " relative=" << sqrt(diff_tot/norm2(rhs,rb[1])) << std::endl;
}
  
  
    
  


  



#if  0
// Check both symm and asymm linops have same unprec op
TEST_F(SymmFixture, CheckUnprecOp)
{
	T x;
	T unprec_symm_x;
	T unprec_asymm_x;

	gaussian(x);
	{

		QDPIO::cout << "Regular Op:" << std::endl;

		(*M_symm).unprecLinOp( unprec_symm_x, x, PLUS);
		(*M_asymm).unprecLinOp( unprec_asymm_x, x, PLUS);

		unprec_asymm_x -= unprec_symm_x;
		Double norm_cb0 = sqrt(norm2(unprec_asymm_x,rb[0]));
		Double norm_cb1 = sqrt(norm2(unprec_asymm_x,rb[1]));

		QDPIO::cout << "CB=0: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb0 <<std::endl;
		QDPIO::cout << "CB=1: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb1 << std::endl;

		ASSERT_LT( toDouble(norm_cb0), 1.0e-13);
		ASSERT_LT( toDouble(norm_cb1), 1.0e-13);

	}

	{
		QDPIO::cout << "Daggered Op:" << std::endl;

		(*M_symm).unprecLinOp( unprec_symm_x, x, MINUS);
		(*M_asymm).unprecLinOp( unprec_asymm_x, x, MINUS);
		unprec_asymm_x -= unprec_symm_x;

		Double norm_cb0 = sqrt(norm2(unprec_asymm_x,rb[0]));
		Double norm_cb1 = sqrt(norm2(unprec_asymm_x,rb[1]));

		QDPIO::cout << "CB=0: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb0 <<std::endl;
		QDPIO::cout << "CB=1: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb1 << std::endl;

		ASSERT_LT( toDouble(norm_cb0), 1.0e-13);
		ASSERT_LT( toDouble(norm_cb1), 1.0e-13);

	}

}

// Check QProp Functionality.
TEST_P(QPropTest, CheckQprop)
{
	LatticeFermion rhs=zero;
	gaussian(rhs);

	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
    Handle<SystemSolver<T>>	qprop_solver = S_symm->qprop(state,inv_param);
	LatticeFermion x = zero;

	(*qprop_solver)(x,rhs);

	// Check residuum
	LatticeFermion Ax=zero;
	(*M_symm).unprecLinOp(Ax,x,PLUS);
	Ax -= rhs;

	Double resid_cb0 = sqrt(norm2(Ax,rb[0]));
	Double resid_cb1 = sqrt(norm2(Ax,rb[1]));
	QDPIO::cout << "Qprop: rsd cb0 = " << resid_cb0 << std::endl;
	QDPIO::cout << "Qprop: rsd cb1 = " << resid_cb1 << std::endl;

	Double resid = sqrt(norm2(Ax));
	Double resid_rel = resid/sqrt(norm2(rhs));
	QDPIO::cout << "QProp Check Back: || r || = " << resid << "  || r ||/||b|| = "
			<< resid_rel << std::endl;

	ASSERT_LT(toDouble(resid_rel), 1.0e-8);

}

INSTANTIATE_TEST_CASE_P(PropSyssolver,
                        QPropTest,
                        ::testing::Values(inv_param_syssolver_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(PropQUDASolver,
				        QPropTest,
						::testing::Values(inv_param_quda_bicgstab_xml,
										  inv_param_quda_multigrid_xml));

#endif

class MdagMInvTestSymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};
class MdagMInvTestAsymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};

TEST_P(MdagMInvTestSymm, CheckMdagMInvSymm)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_symm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;

	(*MdagM_solver)(x,b);

	T tmp, r;

	(*M_symm)(tmp,x,PLUS);
	(*M_symm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);



}

#ifdef BUILD_QUDA
TEST_P(MdagMInvTestSymm, CheckMdagMInvSymmQUDAPredict)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_symm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;
	QUDA4DChronoPredictor chrono(5,DEFAULT);
	(*MdagM_solver)(x,b,chrono);

	T tmp, r;

	(*M_symm)(tmp,x,PLUS);
	(*M_symm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);



}
#endif

TEST_P(MdagMInvTestAsymm, CheckMdagMInvAsymm)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_asymm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;

	(*MdagM_solver)(x,b);

	T tmp, r;

	(*M_asymm)(tmp,x,PLUS);
	(*M_asymm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);



}

#ifdef BUILD_QUDA
TEST_P(MdagMInvTestAsymm, CheckMdagMInvAsymmQUDAPredict)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_asymm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;

	QUDA4DChronoPredictor chrono(5,DEFAULT);

	(*MdagM_solver)(x,b,chrono);

	T tmp, r;

	(*M_asymm)(tmp,x,PLUS);
	(*M_asymm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),2.0e-8);

}
#endif

INSTANTIATE_TEST_CASE_P(MdagMInvSysSolver,
                        MdagMInvTestSymm,
                        ::testing::Values(inv_param_syssolver_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MdagMInvQUDASolver,
                        MdagMInvTestSymm,
						::testing::Values(inv_param_quda_bicgstab_xml,
																  inv_param_quda_multigrid_xml));
#endif

INSTANTIATE_TEST_CASE_P(MdagMInvSysSolver,
                        MdagMInvTestAsymm,
                        ::testing::Values(inv_param_syssolver_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MdagMInvQUDASolver,
                        MdagMInvTestAsymm,
						::testing::Values(inv_param_quda_bicgstab_asymm_xml,
																  inv_param_quda_multigrid_asymm_xml));
#endif

class multiMdagMInvTestSymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};
class multiMdagMInvTestAsymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};

TEST_P(multiMdagMInvTestSymm, checkMultiShift)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");

	Handle<MdagMMultiSystemSolver<T>> multiMdagM = S_symm->mInvMdagM(state,inv_param);

	T rhs;
	gaussian(rhs,rb[1]);

	int n_shift=3;
	multi1d<Real> shifts(n_shift);
	shifts[0] = 0.0001;
	shifts[1] = 0.01;
	shifts[2] = 0.1;

	// Zero the initial guesses
	multi1d<T> solns(n_shift);
	for(int shift=0; shift < n_shift; ++shift) {
		(solns[shift])[rb[1]]=zero;
	}

	//operator() (multi1d< multi1d<T> >& psi, const multi1d<Real>& shifts, const multi1d<T>& chi)
	(*multiMdagM)(solns,shifts,rhs);

	for(int shift = 0; shift < n_shift; ++shift) {
		T r = zero;
		T tmp = zero;
		// r = M^\dag M solns[shift]
		(*M_symm)(tmp,solns[shift],PLUS);
		(*M_symm)(r, tmp, MINUS);

		// r = M^\dag M solns[shift] + shifts[shift] solns[shift]
		//   = (M^\dag M + shifts[shift]) solns[shift]
		r[rb[1]] += shifts[shift]*solns[shift];

		// -residudum
		r[rb[1]] -= rhs;

		Double resid_rel = sqrt(norm2(r,rb[1])/norm2(rhs,rb[1]));
		QDPIO::cout << "shift="<<shift << " || r || / || b ||=" << resid_rel << std::endl;
		ASSERT_LT( toDouble(resid_rel), 1.0e-8);
	}
}

TEST_P(multiMdagMInvTestAsymm, checkMultiShift)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");

	Handle<MdagMMultiSystemSolver<T>> multiMdagM = S_asymm->mInvMdagM(state,inv_param);

	T rhs;
	gaussian(rhs,rb[1]);

	int n_shift=3;
	multi1d<Real> shifts(n_shift);
	shifts[0] = 0.0001;
	shifts[1] = 0.01;
	shifts[2] = 0.1;

	// Zero the initial guesses
	multi1d<T> solns(n_shift);
	for(int shift=0; shift < n_shift; ++shift) {
		(solns[shift])[rb[1]]=zero;
	}

	//operator() (multi1d< multi1d<T> >& psi, const multi1d<Real>& shifts, const multi1d<T>& chi)
	(*multiMdagM)(solns,shifts,rhs);

	for(int shift = 0; shift < n_shift; ++shift) {
		T r = zero;
		T tmp = zero;
		// r = M^\dag M solns[shift]
		(*M_asymm)(tmp,solns[shift],PLUS);
		(*M_asymm)(r, tmp, MINUS);

		// r = M^\dag M solns[shift] + shifts[shift] solns[shift]
		//   = (M^\dag M + shifts[shift]) solns[shift]
		r[rb[1]] += shifts[shift]*solns[shift];

		// -residudum
		r[rb[1]] -= rhs;

		Double resid_rel = sqrt(norm2(r,rb[1])/norm2(rhs,rb[1]));
		QDPIO::cout << "shift="<<shift << " || r || / || b ||=" << resid_rel << std::endl;
		ASSERT_LT( toDouble(resid_rel), 1.0e-8);
	}
}

INSTANTIATE_TEST_CASE_P(MultiShiftSysSolver,
                        multiMdagMInvTestAsymm,
                        ::testing::Values(inv_param_multi_cg_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MultiShiftQUDASolver,
                        multiMdagMInvTestAsymm,
						::testing::Values(inv_param_multi_cg_quda_asymm_xml));
#endif

INSTANTIATE_TEST_CASE_P(MultiShiftSysSolver,
                        multiMdagMInvTestSymm,
                        ::testing::Values(inv_param_multi_cg_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MultiShiftQUDASolver,
                        multiMdagMInvTestSymm,
						::testing::Values(inv_param_multi_cg_quda_xml));
#endif


// Forces
TEST_F(SymmFixture, TestDeriv)
{
	// M_symm = M_oo^{-1} M_asymm
	// X^\dagger d[ M_symm ] Y
	//  = X^\dagger d[ M_oo^{-1} ] M_asymm Y
	//   +X^\dagger M_oo^{-1} d[ M_asymm ] Y
    //
	// = -X^\dagger M_oo^{-1} d[ M_oo ] M^{-1}_oo M_asymm Y
	//   +X^\dagger M_oo^{-1} d[ M_asymm   ] Y
	//
	// =  -Z^\dagger d[ M_oo ] W
	//    +Z^\dagger d[ M_asymm ] Y
	//
	// with W = M_{oo}^{-1} M_asymm Y
	// and  Z = M_oo^{-dagger} X

	P ds_symm;
	P ds_tmp;
	P rhs;


	T X = zero;
	T Y = zero;

	gaussian(X,rb[1]);
	gaussian(Y,rb[1]);

	T tmp = zero;
	T W = zero;
	T Z = zero;

	// W = M^{-1}_oo M_asymm Y
	(*M_asymm)(tmp,Y,PLUS);
	(*M_symm).unprecOddOddInvLinOp(W,tmp,PLUS);

	(*M_symm).unprecOddOddInvLinOp(Z,X,MINUS);

	// The derivative of M_symm
	(*M_symm).deriv(ds_symm,X,Y,PLUS);

	// rhs = Z^\dagger d[ M_asymm ] Y
	(*M_asymm).deriv(rhs,Z,Y,PLUS);

	// rhs -= Z^\dagger d[ M_oo ] W
	(*M_symm).derivUnprecOddOddLinOp(ds_tmp,Z,W,PLUS);
	rhs -= ds_tmp;

	for(int mu=0; mu < Nd; ++mu) {
		rhs[mu] -= ds_symm[mu];
		Double norm_rhs = sqrt(norm2(rhs[mu]));
		Double norm_rhs_per_number = norm_rhs/Double(3*3*2*Layout::vol());
		QDPIO::cout << "mu=" << mu << " || rhs - ds_symm || = " << norm_rhs
				<< "  || rhs - ds_symm || / number =" << norm_rhs_per_number << std::endl;

		ASSERT_LT(toDouble(norm_rhs_per_number), 1.0e-18 );
	}
}

TEST_F(SymmFixture, TestDerivDagger)
{
	// M_symm^\dagger = M_asymm^\dagger M_oo^{-\dagger}

	// X^\dagger d[ M_symm ] Y
	//  = X^\dagger d[ M_asymm^\dagger ] M_oo^{-\dagger} Y
	//   +X^\dagger M_asymm^\dagger d[ M_oo^{-\dagger} ] Y
    //
	// = +X^\dagger d[M_asymm^\dagger ] M_oo^{-dagger} Y
	//   -X^\dagger M_asymm^\dagger ] M_oo^{-dagger} d[ M_oo^{\dagger} ] M_oo^{-dagger} Y
	//
	// =  +X^\dagger d[ M_asymm^\dagger ] W
	//    -Z^\dagger d[ M_oo^\dagger  ] W
	//
	// with W = M_{oo}^{-\dagger}  Y
	// and  Z = M_oo^{-1} M_asymm  X

	P ds_symm;
	P ds_tmp;
	P rhs;


	T X = zero;
	T Y = zero;

	gaussian(X,rb[1]);
	gaussian(Y,rb[1]);

	T tmp = zero;
	T W = zero;
	T Z = zero;

	// W = M^{-dagger}_oo Y
	(*M_symm).unprecOddOddInvLinOp(W,Y,MINUS);

	//Z = M_oo^{-1} M_asymm  X
	(*M_asymm)(tmp,X,PLUS);
	(*M_symm).unprecOddOddInvLinOp(Z,tmp,PLUS);

	// The derivative of M_symm^dagger
	(*M_symm).deriv(ds_symm,X,Y,MINUS);

	// rhs = X^\dagger d[ M_asymm ] W
	(*M_asymm).deriv(rhs,X,W,MINUS);

	// rhs -= Z^\dagger d[ M_oo ] W
	(*M_symm).derivUnprecOddOddLinOp(ds_tmp,Z,W,MINUS);
	rhs -= ds_tmp;

	for(int mu=0; mu < Nd; ++mu) {
		rhs[mu] -= ds_symm[mu];
		Double norm_rhs = sqrt(norm2(rhs[mu]));
		Double norm_rhs_per_number = norm_rhs/Double(3*3*2*Layout::vol());
		QDPIO::cout << "mu=" << mu << " || rhs - ds_symm || = " << norm_rhs
				<< "  || rhs - ds_symm || / number =" << norm_rhs_per_number << std::endl;

		ASSERT_LT(toDouble(norm_rhs_per_number), 1.0e-18 );
	}
}

TEST_F(SymmFixture,TestLogDetUnitGauge)
{
	Q u_unit;
	u_unit.resize(Nd);
	for(int mu=0; mu < Nd; ++mu) {
		u_unit[mu] = 1;
	}

	Handle< FermState<T,P,Q> > unit_state = S_symm->createState(u_unit);
	QDPIO::cout << "Unit Gauge Symm Linop creation" << std::endl;
	Handle<LinOpSymm_T> lop = dynamic_cast<LinOpSymm_T*>(S_symm->linOp(unit_state));

	Double S_e = lop->logDetEvenEvenLinOp();
	Double S_o = lop->logDetOddOddLinOp();
	QDPIO::cout << "Unit gauge: log Det EvenEven = " << S_e << std::endl;
	QDPIO::cout << "Unit gauge: log Det OddOdd = " << S_o << std::endl;

	Double absdiff = fabs(S_e-S_o);
	QDPIO::cout << "Unit gauge: | S_e - S_o | = " << absdiff << std::endl;
	ASSERT_LT( toDouble(absdiff), 1.0e-14);
}

TEST_F(SymmFixture,TestLogDetShiftedGauge)
{
	Q u_shifted;
	u_shifted.resize(Nd);
	for(int mu=0; mu < Nd; ++mu) {
		u_shifted[mu] = shift(u[mu],FORWARD,3);
	}

	Handle< FermState<T,P,Q> > shifted_state = S_symm->createState(u_shifted);

	QDPIO::cout << "Shifted Gauge Symm Linop creation" << std::endl;

	Handle<LinOpSymm_T> shifted_lop = dynamic_cast<LinOpSymm_T*>(S_symm->linOp(shifted_state));

	Double S_e_asymm = M_asymm->logDetEvenEvenLinOp();
	Double S_e = M_symm->logDetEvenEvenLinOp();
	Double S_o = M_symm->logDetOddOddLinOp();

	Double S_e_shifted = shifted_lop->logDetEvenEvenLinOp();
	Double S_o_shifted = shifted_lop->logDetOddOddLinOp();

	QDPIO::cout << "Asymm op log Det EvenEven = " << S_e_asymm << std::endl;
	QDPIO::cout << "Random gauge: log Det EvenEven = " << S_e << std::endl;
	QDPIO::cout << "Random gauge: log Det OddOdd = " << S_o << std::endl;
	QDPIO::cout << "Shifted gauge: log Det EvenEven = " << S_e_shifted << std::endl;
	QDPIO::cout << "Shifted gauge: log Det OddOdd = " << S_o_shifted << std::endl;

	Double diff_e_symm_asymm = fabs(S_e_asymm - S_e);
	Double diff_eo = fabs(S_e - S_o_shifted)/Double(rb[0].numSiteTable());
	Double diff_oe = fabs(S_o - S_e_shifted)/Double(rb[1].numSiteTable());
	QDPIO::cout << "| logDet_e_asymm - logdet_e_symm | = " << diff_e_symm_asymm << std::endl;
	QDPIO::cout << "| logDet_e - logDet_o_shifted |/site = " << diff_eo << std::endl;
	QDPIO::cout << "| logDet_o - logDet_e_shifted |/site = " << diff_oe << std::endl;

	ASSERT_LT( toDouble(diff_e_symm_asymm), 1.0e-14);
	ASSERT_LT( toDouble(diff_eo), 1.0e-13);
	ASSERT_LT( toDouble(diff_oe), 1.0e-13);

}

TEST_F(SymmFixture,TestTwist)
{

	Real twist=Real(0.05);




	LatticeFermion source;
	gaussian(source,rb[1]);
    LatticeFermion t1, t2;

    {
	(*M_tw)(t1,source,PLUS);
	(*M_symm)(t2,source, PLUS);
	t2[ rb[1] ] += twist*(GammaConst<Ns,Ns*Ns-1>()*timesI(source));

	t2[ rb[1] ] -= t1;
	Double normdiff = sqrt(norm2(t2,rb[1]));
	QDPIO::cout << "PLUS : || M(mu) - ( Mdag + i gamma_5 mu ) || = "
			<< normdiff << std::endl;

	ASSERT_LT( toDouble(normdiff), 1.0e-14);
    }

    {
 	(*M_tw)(t1,source,MINUS);
 	(*M_symm)(t2,source, MINUS);
 	t2[ rb[1] ] -= twist*(GammaConst<Ns,Ns*Ns-1>()*timesI(source));

 	t2[ rb[1] ] -= t1;
 	Double normdiff = sqrt(norm2(t2,rb[1]));
 	QDPIO::cout << "MINUS : || M^dag(mu) - ( Mdag - igamma5 mu ) || = "  << normdiff << std::endl;

 	ASSERT_LT( toDouble(normdiff), 1.0e-14);
     }

    {
    	LatticeFermion mdagm,t3;
    	(*M_tw)(t1,source, PLUS);
    	(*M_tw)(t2,t1,MINUS);

    	// M^\dag M
    	(*M_symm)(t1,source,PLUS);
    	(*M_symm)(mdagm,t1,MINUS);

    	// + mu^2 source
    	mdagm[rb[1]] += (twist*twist)*source;

    	// +i mu M^\dag gamma_5
    	t1[rb[1]] = (GammaConst<Ns,Ns*Ns-1>()*timesI(source));
    	(*M_symm)(t3,t1,MINUS);
    	mdagm[rb[1]] += twist*t3;

    	// -i mu gamma_5 M soure
    	(*M_symm)(t1,source,PLUS);
    	t3[rb[1]] = (GammaConst<Ns,Ns*Ns-1>()*timesI(t1));
    	mdagm[rb[1]] -= twist*t3;


    	mdagm[rb[1]] -= t2;
    	Double normdiff = sqrt(norm2(mdagm,rb[1]));
    	QDPIO::cout << "MDAGM : || M^dag(mu)M(mu) - ( Mdag M + mu^2 + imu[ M^dag g_5 - g_5 M ]) || = "  << normdiff << std::endl;
    	ASSERT_LT( toDouble(normdiff), 1.0e-13);
    }
}

class TrLogForceFixture : public SymmFixtureT<::testing::TestWithParam<enum PlusMinus>>{};
TEST_P(TrLogForceFixture,TestShiftedGaugeTrLnForce)
{
	Q u_shifted;
	u_shifted.resize(Nd);
	for(int mu=0; mu < Nd; ++mu) {
		u_shifted[mu] = shift(u[mu],FORWARD,3);
	}

	Handle< FermState<T,P,Q> > periodic_state = S_symm_periodic->createState(u);
	Handle< FermState<T,P,Q> > periodic_shifted_state = S_symm_periodic->createState(u_shifted);

	QDPIO::cout << "Shifted Gauge Symm Linop creation" << std::endl;

	Handle<LinOpAsymm_T> asymm_periodic_op = dynamic_cast<LinOpAsymm_T*>(S_asymm_periodic->linOp(periodic_state));
	Handle<LinOpSymm_T> symm_periodic_op = dynamic_cast<LinOpSymm_T*>(S_symm_periodic->linOp(periodic_state));
	Handle<LinOpSymm_T> shifted_periodic_op = dynamic_cast<LinOpSymm_T*>(S_symm_periodic->linOp(periodic_shifted_state));

	//! Get the force from the EvenEven Trace Log
	//   virtual void derivLogDetEvenEvenLinOp(P& ds_u, enum PlusMinus isign) const = 0;
	P ds_asymm;
	asymm_periodic_op->derivLogDetEvenEvenLinOp(ds_asymm,GetParam());

	P ds_symm;
	symm_periodic_op->derivLogDetEvenEvenLinOp(ds_symm,GetParam());

	P ds_symm_shifted;
	shifted_periodic_op->derivLogDetOddOddLinOp(ds_symm_shifted, GetParam());

	P ds_unshifted; ds_unshifted.resize(Nd);
	P ds_wrong_shifted; ds_wrong_shifted.resize(Nd);
	for(int mu=0; mu < Nd; ++mu ) {
		ds_unshifted[mu] = shift(ds_symm_shifted[mu], BACKWARD,3);
		ds_wrong_shifted[mu] = shift(ds_symm_shifted[mu], FORWARD, 3);
	}
	for(int mu=0; mu < Nd; ++mu ) {
		Double mu_contrib_asymm = norm2(ds_asymm[mu]);

		Double mu_contrib_symm = norm2(ds_symm[mu]);

		Double mu_contrib_shifted = norm2(ds_symm_shifted[mu]);

		QDPIO::cout << "mu=" << mu << "  asymm force norm="<< mu_contrib_asymm << std::endl;
		QDPIO::cout << "mu=" << mu << "  symm force norm ="<< mu_contrib_symm << std::endl;
		QDPIO::cout << "mu=" << mu << "  shifted force norm="<< mu_contrib_shifted << std::endl;


		Double diff_symm_asymm_ee = fabs(mu_contrib_asymm - mu_contrib_symm );
		Double diff_symm_ee_oo = fabs(mu_contrib_symm - mu_contrib_shifted );

		QDPIO::cout << "mu=" << mu << " | F_asymm_ee - F_symm_ee | = " << diff_symm_asymm_ee << std::endl;
		QDPIO::cout << "mu=" << mu << " | F_ee - F_shifted_oo | = " << diff_symm_ee_oo << std::endl;

		ASSERT_LT( toDouble(diff_symm_asymm_ee), 1.0e-15);
		ASSERT_LT( toDouble(diff_symm_ee_oo), 1.0e-15);

		ds_unshifted[mu] -= ds_symm[mu];
		ds_wrong_shifted[mu] -= ds_symm[mu];
		Double diffnorm_unshifted = sqrt(norm2(ds_unshifted[mu]));
		Double diffnorm_per_site = diffnorm_unshifted/Layout::vol();
		QDPIO::cout << "mu=" << mu << " || F_mu - shift(F_shifted, BACKWARD, mu) || =" << diffnorm_unshifted << std::endl;
		QDPIO::cout << "mu=" << mu << " || F_mu - shift(F_shifted, BACKWARD, mu) ||/site =" << diffnorm_per_site << std::endl;

		ASSERT_LT( toDouble(diffnorm_per_site), 1.0e-14 );

		Double diffnorm_wrong_shifted = sqrt(norm2(ds_wrong_shifted[mu]));
		QDPIO::cout << "mu=" << mu << "  DELIBERATELY WRONG: || F_mu - shift(F_shifted, FORWARD, mu) || =" << diffnorm_wrong_shifted << std::endl;
		ASSERT_GT( toDouble(diffnorm_wrong_shifted), 0.5);
		QDPIO::cout << std::endl;
	}

}

INSTANTIATE_TEST_CASE_P(TrLogForces,
		TrLogForceFixture,
		::testing::Values(PLUS,MINUS));


TEST_F(SymmFixture, CheckDerivMultipole)
{
	int n_poles = 5;
	multi1d<T> X(n_poles);
	multi1d<T> Y(n_poles);


	for(int i=0; i < n_poles; ++i) {
		gaussian(X[i],rb[1]);
		gaussian(Y[i],rb[1]);
	}

	multi1d<LatticeColorMatrix> ds_u(Nd);
	multi1d<LatticeColorMatrix> ds_tmp(Nd);
	multi1d<LatticeColorMatrix> ds_mp(Nd);

	for(int dag=0; dag < 2; ++dag ) {
		QDPIO::cout << "Dagger = " << dag << std::endl;

		enum PlusMinus isign = (dag == 0 ) ? PLUS : MINUS;
		// Isign = PLUS
		// Zero sum
		for(int mu=0; mu < Nd; ++mu) {
			ds_u[mu] = zero;
		}

		// Accumulate forces
		StopWatch swatch;
		swatch.start();
		for(int i=0; i < n_poles; ++i) {
			M_symm->deriv(ds_tmp, X[i],Y[i], isign);
			ds_u += ds_tmp;
		}
		swatch.stop();
		QDPIO::cout << "Individual Poles took: " << swatch.getTimeInSeconds() << " sec." << std::endl;

		// do Multipole version
		swatch.reset(); swatch.start();
		M_symm->derivMultipole(ds_mp, X,Y,isign);
		swatch.stop();
		QDPIO::cout << "Multipole took: " << swatch.getTimeInSeconds() << " sec." << std::endl;

		for(int mu=0; mu < Nd; ++mu) {
			ds_mp[mu] -= ds_u[mu];
			Double normdiff = sqrt(norm2(ds_mp[mu]));
			Double normdiff_per_link = normdiff/static_cast<double>(4*Layout::vol());
			QDPIO::cout << "mu="<< mu << " normdiff = " << normdiff
					<< " normdiff per link = " << normdiff_per_link << std::endl;
			ASSERT_LT( toDouble(normdiff_per_link), 1.0e-17);
		}
	}
}

#endif
