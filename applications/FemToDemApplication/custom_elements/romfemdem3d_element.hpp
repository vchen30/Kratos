#if !defined(KRATOS_ROMFEMDEM3D_ELEMENT_H_INCLUDED )
#define  KRATOS_ROMFEMDEM3D_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "custom_elements/femdem3d_element.hpp"

namespace Kratos
{

	class RomFemDem3DElement : public FemDem3DElement    
	{
	public:

		/// Default constructors
		RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry);

		RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		///Copy constructor
		RomFemDem3DElement(RomFemDem3DElement const& rOther);

		/// Destructor.
		virtual ~RomFemDem3DElement();

		/// Assignment operator.
		RomFemDem3DElement& operator=(RomFemDem3DElement const& rOther);


		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


		Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

		RomFemDem3DElement()
		{
		}

        void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);
        void CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    private:












    };

}// Namespace Kratos
#endif // KRATOS_ROMFEMDEM3D_ELEMENT_H_INCLUDED  defined 