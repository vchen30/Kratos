
#if !defined(TEMPLATED_DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  TEMPLATED_DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "includes/serializer.h"
#include "containers/flags.h"

#include "../custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"


namespace Kratos {
    
    class Properties;
    //class TParticleType; // forward declaration of spheric cont particle
    class DEMWall; //forward declaration

    template<class TParticleType>
    class /*__declspec( dllexport )*/ TemplatedDEMDiscontinuumConstitutiveLaw : public Flags {
    public:

        double mKn;
        double mKt;

        KRATOS_CLASS_POINTER_DEFINITION(TemplatedDEMDiscontinuumConstitutiveLaw);

        TemplatedDEMDiscontinuumConstitutiveLaw();

        TemplatedDEMDiscontinuumConstitutiveLaw(const TemplatedDEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& r_process_info);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
        
        virtual std::string GetTypeOfLaw();

        virtual ~TemplatedDEMDiscontinuumConstitutiveLaw();

        virtual TemplatedDEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateContactArea(double radius, double other_radius, double &calculation_area);

        virtual void CalculateElasticConstants(double &kn_el,
                double &kt_el,
                double initial_dist,
                double equiv_young,
                double equiv_poisson,
                double calculation_area,
                TParticleType* element1,
                TParticleType* element2);

        
        virtual void CalculateElasticEnergy(double& normal_elastic_energy,
                                                                double indentation,
                                                                double& cohesive_force,
                                                                TParticleType* element1,
                                                                TParticleType* element2);

        
        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding);

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                TParticleType* element1,
                TParticleType* element2,
                double kn_el,
                double kt_el);

        virtual void InitializeContact(TParticleType * const element1, TParticleType * const element2, const double ini_delta = 0.0);
        virtual void InitializeContactWithFEM(TParticleType* const element, DEMWall* const wall, const double indentation, const double ini_delta = 0.0);
        
        virtual void GetContactStiffness(TParticleType* const element1, TParticleType* const element2, const double ini_delta, double& kn,double& kt);
        
        virtual void CalculateForces(const ProcessInfo& r_process_info,
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
                                                        bool& sliding, double LocalCoordSystem[3][3]);
        
        virtual void CalculateForcesWithFEM( ProcessInfo& r_process_info,
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
                                            bool& sliding);
                
        virtual double CalculateNormalForce(const double indentation);
        virtual double CalculateNormalForce(TParticleType* const element1, TParticleType* const element2, const double indentation, double LocalCoordSystem[3][3]);
        virtual double CalculateNormalForce(TParticleType* const element, DEMWall* const wall, const double indentation);
        virtual double CalculateCohesiveNormalForce(TParticleType * const element1, TParticleType * const element2, const double indentation);
        virtual double CalculateCohesiveNormalForceWithFEM(TParticleType* const element, DEMWall* const wall, const double indentation);
        virtual double LocalPeriod(const int i, TParticleType* element1,TParticleType* element2);


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);

        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

