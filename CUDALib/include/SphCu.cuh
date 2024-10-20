#pragma once
struct Particle {
    double3 Position;
    double3 force;
    double density;
    double3 velocity;
    double3 accel;
    double pressure;
    uint16_t hash;
};
struct SPHSettings
{
    SPHSettings() {}
    SPHSettings(
        double mass, double restDensity, double gasConst, double viscosity,
        double h, double g, double tension)
        : mass(mass)
        , restDensity(restDensity)
        , gasConstant(gasConst)
        , viscosity(viscosity)
        , h(h)
        , g(g)
        , tension(tension) {

        poly6 = 315.0f / (64.0f * 3.14159265359 * pow(h, 9));
        spikyGrad = -45.0f / (3.14159265359 * pow(h, 6));
        spikyLap = 45.0f / (3.14159265359 * pow(h, 6));
        h2 = h * h;
        selfDens = mass * poly6 * pow(h, 6);
        massPoly6Product = mass * poly6;
        sphereScale = h / 2.f;
    };

    double sphereScale;

    double poly6, spikyGrad, spikyLap, gasConstant, mass, h2, selfDens,
        restDensity, viscosity, h, g, tension, massPoly6Product;
};

