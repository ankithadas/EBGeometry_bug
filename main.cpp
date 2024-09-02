#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

// Our include
#include "../../EBGeometry.hpp"

using namespace amrex;

template <class T, class Meta, class BV, size_t K>
class AMReXSDF
{
public:
  AMReXSDF(const std::string a_filename, const bool a_use_bvh)
  {
    if (a_use_bvh) {
      m_sdf = EBGeometry::Parser::readIntoLinearBVH<T, Meta, BV, K>(a_filename);
    }
    else {
      m_sdf = EBGeometry::Parser::readIntoMesh<T, Meta>(a_filename);
    }
  }

  AMReXSDF(const AMReXSDF& a_other)
  {
    this->m_sdf = a_other.m_sdf;
  }

  Real
  operator()(AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
  {
    using Vec3 = EBGeometry::Vec3T<T>;

    return m_sdf->value(Vec3(x, y, z));
  };

  inline Real
  operator()(const RealArray& p) const noexcept
  {
    return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
  }

protected:
  std::shared_ptr<EBGeometry::ImplicitFunction<T>> m_sdf;
};

int
main(int argc, char* argv[])
{

  amrex::Initialize(argc, argv);

  {
    int  max_grid_size   = 16;
    int  which_geom      = 0;
    int  num_coarsen_opt = 0;
    bool use_bvh         = true;

    std::string filename;

    // read parameters
    ParmParse pp;
    pp.query("use_bvh", use_bvh);
    pp.query("max_grid_size", max_grid_size);
    pp.query("which_geom", which_geom);
    pp.query("num_coarsen_opt", num_coarsen_opt);

    Geometry geom;
    {
      RealBox rb;

      rb = RealBox({-20.0, -50.0, 0.}, {332, 238, 160});
      filename = "buildings_fixed.stl"; 

      Array<int, AMREX_SPACEDIM> is_periodic{false, false, false};
      Geometry::Setup(&rb, 0, is_periodic.data());
      Box domain(IntVect(0), IntVect(176-1,144-1,80 - 1));
      geom.define(domain);
    }

    // Create our signed distance function. K is the tree degree while T is the
    // EBGeometry precision.
    constexpr int K = 2;

    using T    = float;
    using Meta = EBGeometry::DCEL::DefaultMetaData;
    using Vec3 = EBGeometry::Vec3T<T>;
    using BV   = EBGeometry::BoundingVolumes::AABBT<T>;

    AMReXSDF<T, Meta, BV, K> sdf(filename, use_bvh);

    auto gshop = EB2::makeShop(sdf);
    EB2::Build(gshop, geom, 0, 0, 1, true, true, num_coarsen_opt);

    // Put some data
    MultiFab mf;
    {
      BoxArray boxArray(geom.Domain());
      boxArray.maxSize(max_grid_size);
      DistributionMapping dm{boxArray};

      std::unique_ptr<EBFArrayBoxFactory> factory =
        amrex::makeEBFabFactory(geom, boxArray, dm, {2, 2, 2}, EBSupport::full);

      mf.define(boxArray, dm, 1, 0, MFInfo(), *factory);
      mf.setVal(1.0);
    }

    EB_WriteSingleLevelPlotfile("plt", mf, {"rho"}, geom, 0.0, 0);
  }

  amrex::Finalize();
}
