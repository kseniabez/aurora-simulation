#pragma once
#include <SFML/Graphics.hpp>
#include <vector>

class FluidSim {
public:
    FluidSim(int width, int height);
    void step(float dt);
    void render(sf::RenderWindow& window, float cellSize);
    float getDensity(int x, int y) const;

private:
    int width, height;
    std::vector<float> density;
};