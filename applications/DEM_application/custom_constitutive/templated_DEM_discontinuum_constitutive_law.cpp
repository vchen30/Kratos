// Last modified by: S. Latorre (CIMNE)
// Date: October 2015

#include "templated_DEM_discontinuum_constitutive_law.h"
#include "custom_elements/templated_spheric_particle.h"
#include "custom_elements/analytic_spheric_particle.h"

namespace Kratos {
    template < class TParticleType >
    TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::TemplatedDEMDiscontinuumConstitutiveLaw() {
        //std::cout << " DEMDiscontinuumConstitutiveLaw constructor..." << std::endl;

    } // Class DEMDiscontinuumConstitutiveLaw
    template < class TParticleType >
    TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::TemplatedDEMDiscontinuumConstitutiveLaw(const TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType> &rReferenceDiscontinuumConstitutiveLaw) {
        //std::cout << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    }

    //    TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::DEMDiscontinuumConstitutiveLaw( const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw){
    //        //std::cout << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    //    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::Initialize(const ProcessInfo& r_process_info) {
    }
//    template < class TParticleType >
//    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
//        //std::cout << "Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
//        pProp->SetValue(TEMPLATED_DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
//    }
    template < class TParticleType >
    std::string TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::GetTypeOfLaw() {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::GetTypeOfLaw) should not be called.","")
        std::string type_of_law = "";
        return type_of_law;
    }
    template < class TParticleType >
    typename TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::Pointer TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::Clone() const {
        typename TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::Pointer p_clone(new TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>(*this));
        return p_clone;
    }
    template < class TParticleType >
    TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::~TemplatedDEMDiscontinuumConstitutiveLaw() {
        //std::cout << "Law destructor..." ; 
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateContactArea(double radius, double other_radius, double &calculation_area) {
        
        KRATOS_TRY
        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = KRATOS_M_PI * rmin*rmin;
        KRATOS_CATCH("")  
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateElasticConstants(double &kn_el,
            double &kt_el,
            double initial_dist,
            double equiv_young,
            double equiv_poisson,
            double calculation_area,
            TParticleType* element1,
            TParticleType* element2) {
        
        KRATOS_TRY 
        double equiv_shear = equiv_young / (2.0 * (1.0 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")  
    }    
    
    template < class TParticleType >
    double TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::LocalPeriod(const int i, TParticleType* element1,
                                               TParticleType* element2) {

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        CalculateContactArea(my_radius, other_radius, calculation_area);

        double radius_sum = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;
        //double initial_delta = element1->GetInitialDelta(i);
        //double initial_dist = radius_sum - initial_delta;

        // calculation of elastic constants for discontinuum - Linear model
        //double kn_el = equiv_young * calculation_area / initial_dist;

        const double modified_radius = equiv_radius * 0.31225; // sqrt(alpha * (2.0 - alpha)) = 0.31225
        double kn = equiv_young * KRATOS_M_PI * modified_radius;        // 2.0 * equiv_young * sqrt_equiv_radius;

        //mDensity = element1_props[PARTICLE_DENSITY];
        //other_density = element2_props[PARTICLE_DENSITY];
        //double m1 = 4/3 * KRATOS_M_PI * my_radius * my_radius * my_radius * mDensity;
        //double m2 = 4/3 * KRATOS_M_PI * other_radius * other_radius * other_radius * other_density;

        const double mRealMass = element1->GetMass();
        const double other_real_mass = element2->GetMass();
        double equiv_mass = (mRealMass*other_real_mass)/(mRealMass+other_real_mass);

        // calculation of damping gamma
        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double friction_coeff = element1->GetProperties()[CONTACT_INTERNAL_FRICC];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);
        const double viscous_damping_coeff     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn);
        double rescaled_damping = viscous_damping_coeff/(2*equiv_mass);
        //double sqr_period = kn / equiv_mass - rescaled_damping*rescaled_damping;
        double sqr_period = sqrt(1+friction_coeff*friction_coeff) * kn / equiv_mass - rescaled_damping*rescaled_damping;
        return sqr_period;
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::InitializeContact(TParticleType* const element1, TParticleType* const element2, const double ini_delta) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::InitializeContact) should not be called.","")
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::InitializeContactWithFEM(TParticleType* const element, DEMWall* const wall, const double indentation, const double ini_delta) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::InitializeContactWithFEM) should not be called.","")
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::GetContactStiffness(TParticleType* const element1, TParticleType* const element2, const double ini_delta, double& kn,double& kt){
      
      InitializeContact(element1, element2, ini_delta);
      kn = mKn;
      kt = mKt;
      
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateForces(const ProcessInfo& r_process_info,
                                                        const double OldLocalContactForce[3],
                                                        double LocalElasticContactForce[3],
                                                        double LocalDeltDisp[3],
                                                        double LocalRelVel[3],            
                                                        double indentation,
                                                        double previous_indentation,
                                                        double ViscoDampingLocalContactForce[3],
                                                        double& cohesive_force,
                                                        TParticleType* element1,
                                                        TParticleType* element2,
                                                        bool& sliding, double LocalCoordSystem[3][3]) {

        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateForces) should not be called.","")
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateElasticEnergy(double& normal_elastic_energy,
                                                                double indentation,
                                                                double& cohesive_force,
                                                                TParticleType* element1,
                                                                TParticleType* element2) {
        
        
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                                                const double OldLocalContactForce[3],
                                                                double LocalElasticContactForce[3],
                                                                double LocalDeltDisp[3],
                                                                double LocalRelVel[3],            
                                                                double indentation,
                                                                double previous_indentation,
                                                                double ViscoDampingLocalContactForce[3],
                                                                double& cohesive_force,
                                                                TParticleType* const element,
                                                                DEMWall* const wall,
                                                                bool& sliding) {
        
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateForcesWithFEM) should not be called.","")
        
    }
    template < class TParticleType >
    double TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateNormalForce(const double indentation) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateNormalForce) should not be called.","")
    }
            
    template < class TParticleType >
    double TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateNormalForce(TParticleType* const element1, TParticleType* const element2, const double indentation,
        double LocalCoordSystem[3][3]) {        
        return CalculateNormalForce(indentation);
    }
    template < class TParticleType >
    double TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateNormalForce(TParticleType* const element, DEMWall* const wall, const double indentation){
        return CalculateNormalForce(indentation);
    }
    template < class TParticleType >
    double TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateCohesiveNormalForce(TParticleType* const element1, TParticleType* const element2, const double indentation) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateCohesiveNormalForce) should not be called.","")
        return 0.0;        
    }
    template < class TParticleType >
    double TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateCohesiveNormalForceWithFEM(TParticleType* const element, DEMWall* const wall, const double indentation){
        KRATOS_THROW_ERROR(std::runtime_error,"This function (TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateCohesiveNormalForceWithFEM) should not be called.","")
        return 0.0;
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateViscoDamping(double LocalRelVel[3],
            double ViscoDampingLocalContactForce[3],
            double indentation,
            double equiv_visco_damp_coeff_normal,
            double equiv_visco_damp_coeff_tangential,
            bool sliding) {

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in opposite direction the visco damping can't overpass the force...

        KRATOS_TRY  

        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];

        if (sliding == false) { //only applied when no sliding to help to the regularized friction law or the spring convergence
            ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        }

        KRATOS_CATCH("")      
    }
    template < class TParticleType >
    void TemplatedDEMDiscontinuumConstitutiveLaw<TParticleType>::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
            double &equiv_visco_damp_coeff_tangential,
            TParticleType* element1,
            TParticleType* element2,
            double kn_el,
            double kt_el) {
        
        KRATOS_TRY 
        double aux_norm_to_tang = 0.0;
        const double my_mass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
        //     const double mDempack_local_damping = element1->GetProperties()[DEMPACK_LOCAL_DAMPING];
        const double mCoefficientOfRestitution = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];

        equiv_visco_damp_coeff_normal = (1-mCoefficientOfRestitution) * 2.0 * sqrt(kn_el / (my_mass + other_real_mass)) * (sqrt(my_mass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // Dempack no ho fa servir...
        KRATOS_CATCH("")  
    }   

    template class TemplatedDEMDiscontinuumConstitutiveLaw<NonTemplatedSphericParticle>;
    template class TemplatedDEMDiscontinuumConstitutiveLaw<AnalyticSphericParticle>;


} // KRATOS
