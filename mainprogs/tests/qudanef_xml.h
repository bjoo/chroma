
/*
 * symm_prec_xml.h
 *
 *  Created on: Oct 18, 2018
 *      Author: bjoo
 */

#ifndef MAINPROGS_TESTS_QUDA_NEF_XML_H_
#define MAINPROGS_TESTS_QUDA_NEF_XML_H_

#include <string>
#include "chroma_config.h"

namespace QUDANefTesting
{
  std::string fermact_xml = \
   "<?xml version='1.0'?>                    \
   <Param>                                   \
    <FermionAction>                          \
      <FermAct>NEF</FermAct>                 \
      <OverMass>1.2</OverMass>               \
      <Mass>0.4</Mass>                       \
      <N5>4</N5>                             \
      <b5>1.2</b5>                           \
      <c5>0.8</c5>                           \
      <AnisoParam>                           \
        <anisoP>false</anisoP>               \
        <t_dir>3</t_dir>                     \
        <xi_0>1</xi_0>                       \
        <nu>1</nu>                           \
      </AnisoParam>                          \
      <FermionBC>                            \
        <FermBC>SIMPLE_FERMBC</FermBC>       \
        <boundary>1 1 1 -1</boundary>        \
      </FermionBC>                           \
    </FermionAction>                         \
  </Param>";

std::string unprec_fermact_xml = \
   "<?xml version='1.0'?>                    \
   <Param>                                   \
    <FermionAction>                          \
      <FermAct>UNPRECONDITIOEND_NEF</FermAct>                 \
      <OverMass>1.2</OverMass>               \
      <Mass>0.4</Mass>                       \
      <N5>4</N5>                             \
      <b5>1.2</b5>                           \
      <c5>0.8</c5>                           \
      <AnisoParam>                           \
        <anisoP>false</anisoP>               \
        <t_dir>3</t_dir>                     \
        <xi_0>1</xi_0>                       \
        <nu>1</nu>                           \
      </AnisoParam>                          \
      <FermionBC>                            \
        <FermBC>SIMPLE_FERMBC</FermBC>       \
        <boundary>1 1 1 -1</boundary>        \
      </FermionBC>                           \
    </FermionAction>                         \
  </Param>";

  
std::string inv_param_quda_cg_xml = \
  "<?xml version='1.0'?>                                                 \
  <Param>					                         \
       <InvertParam>                                                     \
          <invType>QUDA_NEF_INVERTER</invType>                           \
          <NEFParams>                                                    \
            <OverMass>1.2</OverMass>                                     \
            <Mass>0.4</Mass>                                             \
            <N5>4</N5>                                                   \
            <b5>1.2</b5>                                                 \
            <c5>0.8</c5>                                                 \
          </NEFParams>                                                   \
                                                                         \
          <RsdTarget>1.0e-7</RsdTarget>                                  \
          <Delta>1.0e-1</Delta>                                          \
          <MaxIter>1000</MaxIter>                                        \
          <RsdToleranceFactor>100</RsdToleranceFactor>                   \
          <SilentFail>true</SilentFail>                                  \
          <AntiPeriodicT>true</AntiPeriodicT>                            \
          <SolverType>CG</SolverType>                                    \
          <Verbose>false</Verbose>                                       \
          <AsymmetricLinop>true</AsymmetricLinop>                        \
          <CudaReconstruct>RECONS_12</CudaReconstruct>                   \
          <CudaSloppyPrecision>HALF</CudaSloppyPrecision>                \
          <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>       \
          <AxialGaugeFix>false</AxialGaugeFix>                           \
          <AutotuneDslash>true</AutotuneDslash>                          \
        </InvertParam>                                                   \
  </Param>";

std::string inv_param_chroma_cg_xml = \
  "<?xml version='1.0'?>                                                 \
  <Param>					                         \
       <InvertParam>                                                     \
          <invType>CG_INVERTER</invType>                                 \
          <RsdCG>1.0e-7</RsdCG>                                          \
          <MaxCG>1000</MaxCG>                                            \
        </InvertParam>                                                   \
  </Param>";


}// namespace

#endif /* MAINPROGS_TESTS_SYMM_PREC_XML_H_ */
