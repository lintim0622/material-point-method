#include <iostream>
#include "solve.h"

Solve::Solve(mesh_list& meshs)
{
	for (std::unique_ptr<Mesh>& msh : meshs)
    {
		_meshs.push_back(std::move(msh));
	}
	endTime = 0.0;
}

Solve::~Solve()
{

}

void Solve::algorithm(double nowTime)
{
    // for each mesh
    for (std::unique_ptr<Mesh>& msh : _meshs)
    {
        this->calculateParticleInfo(msh);
        this->particleToNode(msh);
        this->nodalSolution(msh);
        this->nodeToParticle(msh, nowTime);
        this->resetNode(msh);
    }
}

void Solve::calculateParticleInfo(std::unique_ptr<Mesh>& msh) 
{
    for (const std::pair<const int, int>& it : msh->pem)
    {
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

void Solve::particleToNode(std::unique_ptr<Mesh>& msh) 
{
    for (const std::pair<const int, int>& it : msh->pem)
    {
        // particle map tp node
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
        Vector2D sInternalForce1{
            ip.dN1[0] * ip.ssp[0] + ip.dN1[1] * ip.ssp[2],
            ip.dN1[0] * ip.ssp[2] + ip.dN1[1] * ip.ssp[1]
        };
        Vector2D sInternalForce2{
            ip.dN2[0] * ip.ssp[0] + ip.dN2[1] * ip.ssp[2],
            ip.dN2[0] * ip.ssp[2] + ip.dN2[1] * ip.ssp[1]
        };
        Vector2D sInternalForce3{
            ip.dN3[0] * ip.ssp[0] + ip.dN3[1] * ip.ssp[2],
            ip.dN3[0] * ip.ssp[2] + ip.dN3[1] * ip.ssp[1]
        };
        Vector2D sInternalForce4{
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
}

void Solve::nodalSolution(std::unique_ptr<Mesh>& msh)
{
    for (Node& node : msh->nodes)
    {
        if (node.mn > 0.0)
        {
            // normalization -> normal vector
            double nNorm{ std::sqrt(std::pow(node.normal[0], 2) + std::pow(node.normal[1], 2)) };
            node.normal *= (1.0 / nNorm);

            // nodal external force
            node.fext = node.bn;

            // nodal total force
            node.ftot = node.fext + node.fint;

            // nodal velocity -> vn
            node.vn = (1.0 / node.mn) * node.pn;

            // nodal acceleration
            node.an = (1.0 / node.mn) * node.ftot;

            // update nodal velocity -> vn+1
            node.vn += DT * node.an;
        }
    }
}

void Solve::nodeToParticle(std::unique_ptr<Mesh>& msh, double nowTime)
{
    for (const std::pair<const int, int>& it : msh->pem)
    {
        // node map to particle
        Particle& ip = msh->particles[it.first];
        Element& ie = msh->elements[it.second];

        const Vector2D& v1{ ie.n1->vn };
        const Vector2D& v2{ ie.n2->vn };
        const Vector2D& v3{ ie.n3->vn };
        const Vector2D& v4{ ie.n4->vn };

        const Vector2D& a1{ ie.n1->an };
        const Vector2D& a2{ ie.n2->an };
        const Vector2D& a3{ ie.n3->an };
        const Vector2D& a4{ ie.n4->an };

        const double& N1{ ip.N1 };
        const double& N2{ ip.N2 };
        const double& N3{ ip.N3 };
        const double& N4{ ip.N4 };

        const Vector2D& dN1{ ip.dN1 };
        const Vector2D& dN2{ ip.dN2 };
        const Vector2D& dN3{ ip.dN3 };
        const Vector2D& dN4{ ip.dN4 };

        ip.xp += DT * (v1 * N1 + v2 * N2 + v3 * N3 + v4 * N4);
        ip.vp += DT * (a1 * N1 + a2 * N2 + a3 * N3 + a4 * N4);

        double depx{ DT * (v1[0] * dN1[0] + v2[0] * dN2[0] + v3[0] * dN3[0] + v4[0] * dN4[0]) };
        double depy{ DT * (v1[1] * dN1[1] + v2[1] * dN2[1] + v3[1] * dN3[1] + v4[1] * dN4[1]) };
        double dexyp = DT * 0.5 * (
            (v1[0] * dN1[1] + v2[0] * dN2[1] + v3[0] * dN3[1] + v4[0] * dN4[1]) +
            (v1[1] * dN1[0] + v2[1] * dN2[0] + v3[1] * dN3[0] + v4[1] * dN4[0])
            );
        ip.ineps[0] = depx;
        ip.ineps[1] = depy;
        ip.ineps[2] = dexyp;

        if (ip.xp[0] > MAXBC || ip.xp[0] < MINBC) {
            std::cout << "t = " << nowTime << std::endl;
            std::cout << "particle " << ip.pid << std::endl;
            std::cout << "xp = [" << std::fixed << std::setprecision(3) 
                      << ip.xp[0] << ", " << ip.xp[1] << "]" << std::endl;
            throw Interpolate_Error("out of boundary for x direction!!!");
        }
        if (ip.xp[1] > MAXBC || ip.xp[1] < MINBC) {
            std::cout << "t = " << nowTime << std::endl;
            std::cout << "particle " << ip.pid << std::endl;
            std::cout << "xp = [" << std::fixed << std::setprecision(3) 
                      << ip.xp[0] << ", " << ip.xp[1] << "]" << std::endl;
            throw Interpolate_Error("out of boundary for y direction!!!");
        }
    }
}

void Solve::updateParticles(std::unique_ptr<Mesh>& msh)
{
    for (const std::pair<const int, int>& it : msh->pem) 
    {
        Particle& ip = msh->particles[it.first];
        const Material& material = *(msh->material);

        double insp[3] {
            material.E1 * ip.ineps[0] + material.E2 * ip.ineps[1], 
            material.E2 * ip.ineps[0] + material.E1 * ip.ineps[1], 
            material.G * ip.ineps[2]
        };
        ip.ep[0] += ip.ineps[0];
        ip.ep[1] += ip.ineps[1];
        ip.ep[2] += ip.ineps[2];

        ip.ssp[0] += insp[0] / material.rho;
        ip.ssp[1] += insp[1] / material.rho;
        ip.ssp[2] += insp[2] / material.rho;
    }

    msh->clearMapId();
    msh->createElementParticleMap();
}

void Solve::resetNode(std::unique_ptr<Mesh>& msh)
{
    for (Node& node : msh->nodes)
    {
        if (node.mn > 0.0)
        {
            node.mn = 0.0;
            node.vn.setZero();
            node.an.setZero();
            node.pn.setZero();
            node.fint.setZero();
            node.fext.setZero();
            node.ftot.setZero();
            node.bn.setZero();
            node.normal.setZero();
        }
    } 
}