#include <iostream>
#include "solve.h"

Solve::Solve(mesh_list& meshs) {
	for (std::unique_ptr<Mesh>& msh : meshs) {
		_meshs.push_back(std::move(msh));
	}
	endTime = 0.0;
}

Solve::~Solve() {

}

void Solve::calculateParticleInfo() {

    // Calculate particle mass and momentum
    for (std::unique_ptr<Mesh>& msh : _meshs) {

        for (const auto& it : msh->pem) { 

            // `it` is of type `const std::pair<const int, int>&`
            Particle& ip = msh->particles[it.first];
            Element& ie = msh->elements[it.second];
            const Material& material = *(msh->material);

            // Calculate particle mass and momentum
            ip.calculateMass(material);
            ip.calculateMomentum();

            // n1 position
            double& xo = ie.n1->xn[0];
            double& yo = ie.n1->xn[1];

            // Particle position
            double& x = ip.xp[0];
            double& y = ip.xp[1];

            // Natural coordinates
            double xi = 2.0 * (x - xo) / UNITGRID - 1.0;
            double eta = 2.0 * (y - yo) / UNITGRID - 1.0;

            // Assign shape functions
            ip.N1 = shapeN1(xi, eta);
            ip.N2 = shapeN2(xi, eta);
            ip.N3 = shapeN3(xi, eta);
            ip.N4 = shapeN4(xi, eta);
            ip.dN1 = { shapedN1x(eta), shapedN1y(xi) };
            ip.dN2 = { shapedN2x(eta), shapedN2y(xi) };
            ip.dN3 = { shapedN3x(eta), shapedN3y(xi) };
            ip.dN4 = { shapedN4x(eta), shapedN4y(xi) };
        }
    }
}

void Solve::particleToNode() {

    // for each mesh
    for (std::unique_ptr<Mesh>& msh : _meshs) {

        // particle map tp node
        for (const auto& it : msh->pem) {
            Particle& ip = msh->particles[it.first];
            Element& ie = msh->elements[it.second];

            // nodal mass-> mn
            ie.n1->mn += ip.mp * ip.N1;
            ie.n2->mn += ip.mp * ip.N2;
            ie.n3->mn += ip.mp * ip.N3;
            ie.n4->mn += ip.mp * ip.N4;

            // nodal momentum -> pn
            ie.n1->pn += ip.Pp * ip.N1;
            ie.n2->pn += ip.Pp * ip.N2;
            ie.n3->pn += ip.Pp * ip.N3;
            ie.n4->pn += ip.Pp * ip.N4;

            // nodal body force
            ie.n1->bn += ip.mp * ip.bp * ip.N1;
            ie.n2->bn += ip.mp * ip.bp * ip.N2;
            ie.n3->bn += ip.mp * ip.bp * ip.N3;
            ie.n4->bn += ip.mp * ip.bp * ip.N4;
            
            // nodal internal force
            Vector2D sInternalForce1 {
                ip.dN1[0] * ip.ssp[0] + ip.dN1[1] * ip.ssp[2],
                ip.dN1[0] * ip.ssp[2] + ip.dN1[1] * ip.ssp[1]
            };
            Vector2D sInternalForce2 {
                ip.dN2[0] * ip.ssp[0] + ip.dN2[1] * ip.ssp[2],
                ip.dN2[0] * ip.ssp[2] + ip.dN2[1] * ip.ssp[1]
            };
            Vector2D sInternalForce3 {
                ip.dN3[0] * ip.ssp[0] + ip.dN3[1] * ip.ssp[2],
                ip.dN3[0] * ip.ssp[2] + ip.dN3[1] * ip.ssp[1]
            };
            Vector2D sInternalForce4 {
                ip.dN4[0] * ip.ssp[0] + ip.dN4[1] * ip.ssp[2],
                ip.dN4[0] * ip.ssp[2] + ip.dN4[1] * ip.ssp[1]
            };
            ie.n1->fint += ip.mp * sInternalForce1;
            ie.n2->fint += ip.mp * sInternalForce2;
            ie.n3->fint += ip.mp * sInternalForce3;
            ie.n4->fint += ip.mp * sInternalForce4;

            // nodal normal vector
            ie.n1->normal += ip.mp * ip.dN1;
            ie.n2->normal += ip.mp * ip.dN2;
            ie.n3->normal += ip.mp * ip.dN3;
            ie.n4->normal += ip.mp * ip.dN4;
        }

        // run all node
        for (Node& node : msh->nodes) {

            if (node.mn > 0.0) {

                // nodal external force
                node.fext = node.bn;

                // normalization
                double nNorm{ std::sqrt(std::pow(node.normal[0], 2) + std::pow(node.normal[1], 2)) };
                node.normal *= (1.0 / nNorm);
            }
        }
    }
}

void Solve::nodalSolution() {

}

void Solve::nodeToParticle() {

}

void Solve::updateParticles() {

}

void Solve::resetNode() {

}