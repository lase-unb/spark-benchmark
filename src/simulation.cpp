#include "simulation.h"

#include <spark/collisions/mcc.h>
#include <spark/constants/constants.h>
#include <spark/core/matrix.h>
#include <spark/em/electric_field.h>
#include <spark/em/poisson.h>
#include <spark/interpolate/field.h>
#include <spark/interpolate/weight.h>
#include <spark/particle/boundary.h>
#include <spark/particle/pusher.h>
#include <spark/random/random.h>
#include <spark/spatial/grid.h>

#include "reactions.h"

#include <fstream>
#include <filesystem>

namespace {
    auto maxwellian_emitter(double t, double lx, double ly, double m) {
        return [t, lx, ly, m](spark::core::Vec<3>& v, spark::core::Vec<2>& x) {
            x.x = lx * spark::random::uniform();
            x.y = ly * spark::random::uniform();
            double vth = std::sqrt(spark::constants::kb * t / m);
            v = {spark::random::normal(0.0, vth), spark::random::normal(0.0, vth),
                 spark::random::normal(0.0, vth)};
        };
    }
} // namespace

namespace spark {
    Simulation::Simulation(const Parameters& parameters, const std::string& data_path)
        : parameters_(parameters), data_path_(data_path), state_(StateInterface(*this)) {}

    void Simulation::run() {
        set_initial_conditions();    
        spatial::TUniformGrid<core::TVec<double, 2>, 2> electric_field_({parameters_.lx, parameters_.ly},{parameters_.nx, parameters_.ny});
        
        auto electron_collisions = load_electron_collisions();
        auto ion_collisions = load_ion_collisions();

        em::StructPoissonSolver2D::DomainProp domain_prop;
        domain_prop.extents = {static_cast<int>(parameters_.nx), static_cast<int>(parameters_.ny)};
        domain_prop.dx = {parameters_.dx, parameters_.dy};

        std::vector<em::StructPoissonSolver2D::Region> region;

        region.push_back({em::CellType::BoundaryDirichlet,
                          {0, 0},
                          {static_cast<int>(parameters_.nx - 1), 0},
                          []()
                          { return 0.0; }});

        region.push_back({em::CellType::BoundaryDirichlet,
                          {0, static_cast<int>(parameters_.ny - 1)},
                          {static_cast<int>(parameters_.nx - 1), static_cast<int>(parameters_.ny - 1)},
                          []()
                          { return 0.0; }});

        region.push_back({em::CellType::BoundaryNeumann,
                          {0, 0},
                          {0, static_cast<int>(parameters_.ny - 1)},
                          []()
                          { return 0.0; }});

        region.push_back({em::CellType::BoundaryNeumann,
                          {static_cast<int>(parameters_.nx - 1), 0},
                          {static_cast<int>(parameters_.nx - 1), static_cast<int>(parameters_.ny - 1)},
                          []()
                          { return 0.0; }});


        auto poisson_solver = em::StructPoissonSolver2D(domain_prop, region);
        events().notify(Event::Start, state_);

        for (step = 0; step < parameters_.n_steps; ++step) {
            spark::interpolate::weight_to_grid(electrons_, electron_density_);
            spark::interpolate::weight_to_grid(ions_, ion_density_);
            
            spark::em::charge_density<2>(parameters_.particle_weight, ion_density_,
                                      electron_density_, rho_field_);


            poisson_solver.solve(phi_field_.data(), rho_field_.data());

            spark::em::electric_field<2>(phi_field_, electric_field_.data());

            spark::interpolate::field_at_particles(electric_field_, electrons_, electron_field);
            spark::interpolate::field_at_particles(electric_field_, ions_, ion_field);

            spark::particle::move_particles(electrons_, electron_field, parameters_.dt);
            spark::particle::move_particles(ions_, ion_field, parameters_.dt);

            tiled_boundary_.apply(&electrons_);
            tiled_boundary_.apply(&ions_);

            electron_collisions.react_all();
            ion_collisions.react_all();

            events().notify(Event::Step, state_);
    } 
            events().notify(Event::End, state_);
}

Events<Simulation::Event, Simulation::EventAction>& Simulation::events() {
    return events_;
}

    void Simulation::set_initial_conditions() {
        electrons_ = spark::particle::ChargedSpecies<2, 3>(-spark::constants::e, spark::constants::m_e);
        electrons_.add(parameters_.n_initial,
                       maxwellian_emitter(parameters_.te, parameters_.lx, parameters_.ly,
                                          spark::constants::m_e));

        ions_ = spark::particle::ChargedSpecies<2, 3>(spark::constants::e, parameters_.m_he);
        ions_.add(parameters_.n_initial,
                  maxwellian_emitter(parameters_.ti, parameters_.lx, parameters_.ly, parameters_.m_he));

        electron_density_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                           {parameters_.nx, parameters_.ny});
        ion_density_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                      {parameters_.nx, parameters_.ny});
        rho_field_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                    {parameters_.nx, parameters_.ny});
        phi_field_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                    {parameters_.nx, parameters_.ny});        

        electron_field = spark::core::TMatrix<spark::core::TVec<double, 2>, 1>(spark::core::ULongVec<1>{parameters_.nx * parameters_.ny});
        ion_field = spark::core::TMatrix<spark::core::TVec<double, 2>, 1>(spark::core::ULongVec<1>{parameters_.nx * parameters_.ny});
        
        std::vector<spark::particle::TiledBoundary> boundaries = {
        {{0, 0}, {static_cast<int>(parameters_.nx - 1), 0}, spark::particle::BoundaryType::Absorbing},
        {{0, static_cast<int>(parameters_.ny - 1)}, {static_cast<int>(parameters_.nx - 1), static_cast<int>(parameters_.ny - 1)}, spark::particle::BoundaryType::Absorbing},
        {{0, 0}, {0, static_cast<int>(parameters_.ny - 1)}, spark::particle::BoundaryType::Absorbing},
        {{static_cast<int>(parameters_.nx - 1), 0}, {static_cast<int>(parameters_.nx - 1), static_cast<int>(parameters_.ny - 1)}, spark::particle::BoundaryType::Absorbing}
    };
        tiled_boundary_ = spark::particle::TiledBoundary2D(electric_field_.prop(), boundaries, parameters_.dt);
    }

    spark::collisions::MCCReactionSet<2, 3> Simulation::load_electron_collisions() {
        auto electron_reactions = reactions::load_electron_reactions(data_path_, parameters_, ions_);
        spark::collisions::ReactionConfig<2, 3> electron_reaction_config{
            parameters_.dt, parameters_.dx,
            std::make_unique<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ng, parameters_.tg),
            std::move(electron_reactions), spark::collisions::RelativeDynamics::FastProjectile};

        return spark::collisions::MCCReactionSet(electrons_, std::move(electron_reaction_config));
    }

    spark::collisions::MCCReactionSet<2, 3> Simulation::load_ion_collisions() {
        auto ion_reactions = reactions::load_ion_reactions(data_path_, parameters_);
        spark::collisions::ReactionConfig<2, 3> ion_reaction_config{
            parameters_.dt, parameters_.dx,
            std::make_unique<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ng, parameters_.tg),
            std::move(ion_reactions), spark::collisions::RelativeDynamics::SlowProjectile};

        return spark::collisions::MCCReactionSet(ions_, std::move(ion_reaction_config));
    }
} // namespace spark
