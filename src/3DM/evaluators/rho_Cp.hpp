//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef _3DM_RHO_CP_HPP
#define _3DM_RHO_CP_HPP

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
  ///
  /// This evaluator computes the specific heat times density
  ///
  template<typename EvalT, typename Traits>
  class rho_Cp : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {

  public:

    rho_Cp(Teuchos::ParameterList& p,
	   const Teuchos::RCP<Albany::Layouts>& dl);

    void 
    postRegistrationSetup(typename Traits::SetupData d,
			  PHX::FieldManager<Traits>& vm);

    void 
    evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename EvalT::ScalarT ScalarT;
    typedef typename EvalT::MeshScalarT MeshScalarT;

    ScalarT constant_value_;

    PHX::MDField<MeshScalarT,Cell,QuadPoint,Dim> coord_;
    PHX::MDField<ScalarT,Cell,QuadPoint> rho_cp_;

    unsigned int num_qps_;
    unsigned int num_dims_;
    unsigned int num_nodes_;
    unsigned int workset_size_;
  
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidrho_CpParameters() const;

  };
}

#endif
