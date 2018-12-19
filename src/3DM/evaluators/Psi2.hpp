//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef _3DM_PSI2_HPP
#define _3DM_PSI2_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Sacado_ParameterAccessor.hpp"
#include "Teuchos_Array.hpp"
#include "Albany_Layouts.hpp"

namespace _3DM {
  template<typename EvalT, typename Traits>
  class Psi2 : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {

  public:

    Psi2(Teuchos::ParameterList& p,
	 const Teuchos::RCP<Albany::Layouts>& dl);

    void 
    postRegistrationSetup(typename Traits::SetupData d,
			  PHX::FieldManager<Traits>& vm);

    void 
    evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename EvalT::ScalarT ScalarT;
    typedef typename EvalT::MeshScalarT MeshScalarT;

    double constant_value_;

    PHX::MDField<ScalarT,Cell,QuadPoint> phi2_;
    PHX::MDField<ScalarT,Cell,QuadPoint> psi2_;


    unsigned int num_qps_;
    unsigned int num_dims_;
    unsigned int num_nodes_;
    unsigned int workset_size_;

    bool enable_transient_;
    std::string psi2_Name_;
    std::string phi2_Name_;
  
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidPsi2Parameters() const; 

  };
}

#endif
