#include <SFML/Graphics.hpp>
#include "fluid_sim.hpp"

int main() {
    const int windowWidth = 800;
    const int windowHeight = 800;

    const int gridWidth = 100;
    const int gridHeight = 100;
    const float cellSize = windowWidth / static_cast<float>(gridWidth);

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Fluid Sim Starter");

    FluidSim sim(gridWidth, gridHeight);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // simulation step (placeholder)
        sim.step(0.016f); // ~60 FPS

        window.clear();
        sim.render(window, cellSize);
        window.display();
    }

    return 0;
}