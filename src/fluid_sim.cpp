#include "fluid_sim.hpp"

FluidSim::FluidSim(int w, int h)
    : width(w), height(h), density(w* h, 0.0f) {
    // Just add a little blob in the center
    int cx = width / 2;
    int cy = height / 2;
    density[cy * width + cx] = 1.0f;
}

void FluidSim::step(float dt) {
    // Just fade the density a little for demo
    for (float& d : density) {
        d *= 0.99f;
    }
}

float FluidSim::getDensity(int x, int y) const {
    return density[y * width + x];
}

void FluidSim::render(sf::RenderWindow& window, float cellSize) {
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float d = getDensity(x, y);
            if (d > 0.01f) {
                sf::RectangleShape rect(sf::Vector2f(cellSize, cellSize));
                rect.setPosition(x * cellSize, y * cellSize);
                rect.setFillColor(sf::Color(0, static_cast<sf::Uint8>(255 * d), 255));
                window.draw(rect);
            }
        }
    }
}
