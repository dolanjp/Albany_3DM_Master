//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//
#include <fstream>
#include "Sacado_ParameterRegistration.hpp"
#include "Albany_Utils.hpp"

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

namespace _3DM {

  //**********************************************************************
  template<typename EvalT, typename Traits>
  ThCond<EvalT, Traits>::
  ThCond(Teuchos::ParameterList& p,
	 const Teuchos::RCP<Albany::Layouts>& dl) :
    coord_      (p.get<std::string>("Coordinate Name"), dl->qp_vector),
    T_          (p.get<std::string>("Temperature Name"), dl->qp_scalar),
    k_          (p.get<std::string>("Thermal Conductivity Name"), dl->qp_scalar),
    phi1_       (p.get<std::string>("Phi1 Name"), dl->qp_scalar),
	psi1_		(p.get<std::string>("Psi1 Name"), dl->qp_scalar),
    psi2_       (p.get<std::string>("Psi2 Name"), dl->qp_scalar)

  {

    this->addDependentField(coord_);
    this->addDependentField(T_);
    this->addDependentField(phi1_);
    this->addDependentField(psi1_);	
    this->addDependentField(psi2_);
	
    this->addEvaluatedField(k_);
 
    Teuchos::RCP<PHX::DataLayout> scalar_dl = dl->qp_scalar;
    std::vector<PHX::Device::size_type> dims;
    scalar_dl->dimensions(dims);
    workset_size_ = dims[0];
    num_qps_      = dims[1];

	
    Teuchos::ParameterList* cond_list =
       p.get<Teuchos::ParameterList*>("Parameter List");

    Teuchos::RCP<const Teuchos::ParameterList> reflist =
       this->getValidThCondParameters();

    cond_list->validateParameters(*reflist, 0,
       Teuchos::VALIDATE_USED_ENABLED, Teuchos::VALIDATE_DEFAULTS_DISABLED);
	   
	//Parameters for the functional form of the thermal conductivity
    a = cond_list->get("a", 1.0);
    b = cond_list->get("b", 1.0);
    c = cond_list->get("c", 1.0);
    d = cond_list->get("d", 1.0);
    e =  cond_list->get("e", 1.0);

    Kl_ =  cond_list->get("Thermal Conductivity Liquid", 1.0);
    Kv_ =  cond_list->get("Thermal Conductivity Vapor", 1.0);
	
	Teuchos::ParameterList* input_list =
		p.get<Teuchos::ParameterList*>("Input List");  
		
	sim_type = input_list->get<std::string>("Simulation Type");
		
	if (sim_type == "SLM Additive"){
		Teuchos::ParameterList* cond_list = p.get<Teuchos::ParameterList*>("Powder Parameter List");
		//Assume constant thermal conductivity in the powder
		Kp_ = cond_list->get("a", 1.0);
	}	
	else{
		Kp_ = 0;
	}
	
    this->setName("ThCond"+PHX::typeAsString<EvalT>());
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  void ThCond<EvalT, Traits>::
  postRegistrationSetup(typename Traits::SetupData d,
			PHX::FieldManager<Traits>& fm)
  {
    this->utils.setFieldData(coord_,fm);
    this->utils.setFieldData(T_,fm);
    this->utils.setFieldData(phi1_,fm);
    this->utils.setFieldData(psi1_,fm);	
    this->utils.setFieldData(psi2_,fm);
    this->utils.setFieldData(k_,fm);
  }

  //**********************************************************************

  template<typename EvalT, typename Traits>
  void ThCond<EvalT, Traits>::
  evaluateFields(typename Traits::EvalData workset)
  {

    // thermal conductivity
    for (std::size_t cell = 0; cell < workset.numCells; ++cell) {
        for (std::size_t qp = 0; qp < num_qps_; ++qp){
            Kd_ = (a + b*T_(cell, qp) + c*T_(cell, qp)*T_(cell, qp) + d/T_(cell, qp) + e/(T_(cell, qp)*T_(cell, qp)));
			Ks_ = (1 - psi_1(cell,qp))*Kp_ + psi_1(cell,qp)*Kd_
            k_(cell, qp) = (Ks_*(1.0 - phi1_(cell, qp)) + Kl_*phi1_(cell, qp))*(1.0 - psi2_(cell, qp)) + Kv_*psi2_(cell, qp);
			//std::cout<<"a = "<< a <<", b = "<< b <<"c = "<< c <<", d = "<< d <<", e ="<< e <<", K = "<< k_(cell, qp) <<std::endl;
	  }
    }
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  Teuchos::RCP<const Teuchos::ParameterList>
  ThCond<EvalT, Traits>::
  getValidThCondParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> valid_pl =
    rcp(new Teuchos::ParameterList("Valid Thermal Conductivity Params"));;

    valid_pl->set<double>("a", 1.0);
    valid_pl->set<double>("b", 1.0);
    valid_pl->set<double>("c", 1.0);
    valid_pl->set<double>("d", 1.0);
    valid_pl->set<double>("e", 1.0);
	
    valid_pl->set<double>("Thermal Conductivity Liquid", 1.0);
    valid_pl->set<double>("Thermal Conductivity Vapor", 1.0);

    return valid_pl;
  }

  //**********************************************************************
}
