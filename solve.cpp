#include <iostream>
#include "solve.h"

// ****************************    SOLVE    ***************************************
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

void Solve::algorithm(double nowTime, std::vector<Boundary>& bcArray,
                      const std::function<double(double)>& decayFunction)
{
    // for each mesh
    for (std::unique_ptr<Mesh>& msh : _meshs)
    {
        this->calculateParticleInfo(msh);
        this->particleToNode(msh);
        this->nodalSolution(msh);
        this->frameBoundary(msh, bcArray, decayFunction);
        this->nodeToParticle(msh, nowTime);
        this->updateParticles(msh);
    }
}

void Solve::calculateParticleInfo(std::unique_ptr<Mesh>& msh) 
{
    for (const std::pair<const int, int>& it : msh->pem)
    {
        // `it` is of type `const std::pair<const int, int>&`
        Particle& ip = msh->particles[it.first];
        Element& ie = msh->elements[it.second];
        const Material& material = msh->material;

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
        ie.n1->fint += (-1.0) * ip.mp * sInternalForce1;
        ie.n2->fint += (-1.0) * ip.mp * sInternalForce2;
        ie.n3->fint += (-1.0) * ip.mp * sInternalForce3;
        ie.n4->fint += (-1.0) * ip.mp * sInternalForce4;

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
            node.normal /= nNorm;

            // nodal external force
            node.fext = node.bn;

            // nodal total force
            node.ftot = node.fext + node.fint;

            // nodal velocity -> vn
            node.vn = node.pn / node.mn;

            // nodal acceleration
            node.an = node.ftot / node.mn;

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
        double depxy = DT * 0.5 * (
            (v1[0] * dN1[1] + v2[0] * dN2[1] + v3[0] * dN3[1] + v4[0] * dN4[1]) +
            (v1[1] * dN1[0] + v2[1] * dN2[0] + v3[1] * dN3[0] + v4[1] * dN4[0])
            );
        ip.ineps[0] = depx;
        ip.ineps[1] = depy;
        ip.ineps[2] = depxy;

        if (ip.xp[0] > MAXBC || ip.xp[0] < MINBC) 
        {
            std::cout << "t = " << nowTime << std::endl;
            std::cout << "particle " << ip.pid << std::endl;
            std::cout << "xp = [" << std::fixed << std::setprecision(6) 
                      << ip.xp[0] << ", " << ip.xp[1] << "]" << std::endl;
            throw Interpolate_Error("out of boundary for x direction!!!");
        }

        if (ip.xp[1] > MAXBC || ip.xp[1] < MINBC) 
        {
            std::cout << "t = " << nowTime << std::endl;
            std::cout << "particle " << ip.pid << std::endl;
            std::cout << "xp = [" << std::fixed << std::setprecision(6) 
                      << ip.xp[0] << ", " << ip.xp[1] << "]" << std::endl;
            throw Interpolate_Error("out of boundary for y direction!!!");
        }
    }
}

static double norm(const Vector2D& vector)
{
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);
}

void Solve::frameBoundary(std::unique_ptr<Mesh>& msh, 
                          std::vector<Boundary>& bcArray, 
                          const std::function<double(double)>& decayFunction)
{
    for (Node& node : msh->nodes)
    {
        if (node.mn > 0.0)
        {
            double& m_ik{ node.mn };
            Vector2D& x_ik{ node.xn };
            Vector2D& vtr_iL{ node.vn };
            Vector2D ptr_iL{ vtr_iL * m_ik };

            for (Boundary& bc : bcArray)
            {
                // projection of a point in global coordinate
                Vector2D x_ikl{ bc.transformGlobal(x_ik) };

                // projection of a point in local coordinate
                Vector2D x_bc{ bc.transformLocal(bc.newtonMethod(x_ikl)) };
                
                // the distance parallel to the normal vector
                double r{ (x_ik - x_bc).dot(bc.nbc()) };

                /* the length of the grid size in the normal direction -> UNITGRID
                   nodal momentum in the normal direction => scale */
                double ptr_iLn{ ptr_iL.dot(bc.nbc()) };

                // influence coefficient->q
                double q{ layer(r, decayFunction) };

                if (bc.form() == "slip")
                {
                    if (ptr_iLn < 0.0 && r < UNITGRID)
                    {
                        // normal force -> (-node.f_int - node.f_ext + node.mass * (-v_ik) / dt).dot(bc.nbc)
                        double fN{ (-ptr_iL / DT).dot(bc.nbc()) };
                        Vector2D f_in{ fN * bc.nbc() };

                        // friction force
                        Vector2D f_iLt{ -(ptr_iL - ptr_iLn * bc.nbc()) / DT };
                        Vector2D f_if{};
                        double norm_f_iLt{ norm(f_iLt) };
                        if (norm_f_iLt != 0.0)
                        {
                            Vector2D et{ f_iLt / norm_f_iLt };
                            f_if = std::min(bc.mu() * fN, norm_f_iLt) * et;
                        }

                        // boundary force
                        Vector2D f_ibc{ q * (f_in + f_if) };
                        Vector2D a_ibc{ f_ibc / m_ik };

                        // update MPM node acceleration && velocity
                        node.fbc += f_ibc;
                        node.an += a_ibc;
                        node.vn += DT * a_ibc;
                    }
                }
                
                else if (bc.form() == "sticky")
                {
                    // boundary force -> q*(-node.f_int - node.f_ext + node.mass*(-v_ik)/dt)
                    Vector2D f_ibc{ q * (-ptr_iL / DT) };
                    Vector2D a_ibc{ f_ibc / m_ik };

                    // update MPM node acceleration && velocity
                    node.fbc += f_ibc;
                    node.an += a_ibc;
                    node.vn += DT * a_ibc;
                }
            }
        }
    }
}

void Solve::updateParticles(std::unique_ptr<Mesh>& msh)
{
    for (const std::pair<const int, int>& it : msh->pem) 
    {
        Particle& ip = msh->particles[it.first];
        const Material& material = msh->material;

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

void Solve::resetNode()
{
    for (std::unique_ptr<Mesh>& msh : _meshs)
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
                node.fbc.setZero();
                node.bn.setZero();
                node.normal.setZero();
            }
        } 
    }
}

void Solve::data_output(const std::string& pfile_name, const std::string& nfile_name, bool append) const {
    std::ios_base::openmode mode = std::ios_base::out;
    if (append) {
        mode |= std::ios_base::app;
    }

    // Particle data output
    std::ofstream pfile(pfile_name, mode);
    if (!pfile.is_open()) {
        std::cerr << "Failed to open particle file: " << pfile_name << std::endl;
        return;
    }

    if (!append) {
        // Assuming particle_info_name contains the headers for the particle file
        pfile << std::setw(8) << "PID"
            << std::setw(14) << "Mass"
            << std::setw(14) << "PosX"
            << std::setw(14) << "PosY"
            << std::setw(14) << "VelX"
            << std::setw(14) << "VelY"
            << std::setw(14) << "StressXX"
            << std::setw(14) << "StressYY"
            << std::setw(14) << "StressXY"
            << std::setw(14) << "StrainXX"
            << std::setw(14) << "StrainYY"
            << std::setw(14) << "StrainXY" << std::endl;
    }

    for (const std::unique_ptr<Mesh>& msh : _meshs) {
        for (const Particle& particle : msh->particles) {
            pfile << std::setw(8) << particle.pid
                << std::setw(14) << std::scientific << particle.mp
                << std::setw(14) << particle.xp[0]
                << std::setw(14) << particle.xp[1]
                << std::setw(14) << particle.vp[0]
                << std::setw(14) << particle.vp[1]
                << std::setw(14) << particle.ssp[0]
                << std::setw(14) << particle.ssp[1]
                << std::setw(14) << particle.ssp[2]
                << std::setw(14) << particle.ep[0]
                << std::setw(14) << particle.ep[1]
                << std::setw(14) << particle.ep[2] << std::endl;
        }
    }
    pfile.close();
    // std::cout << pfile_name << " has been output" << std::endl;

    // Node data output
    std::ofstream nfile(nfile_name, mode);
    if (!nfile.is_open()) {
        std::cerr << "Failed to open node file: " << nfile_name << std::endl;
        return;
    }

    if (!append) {
        // Assuming node_info_name contains the headers for the node file
        nfile << std::setw(8) << "NID"
            << std::setw(14) << "Mass"
            << std::setw(14) << "PosX"
            << std::setw(14) << "PosY"
            << std::setw(14) << "VelX"
            << std::setw(14) << "VelY"
            << std::setw(14) << "AccX"
            << std::setw(14) << "AccY"
            << std::setw(14) << "FIntX"
            << std::setw(14) << "FIntY"
            << std::setw(14) << "FExtX"
            << std::setw(14) << "FExtY" 
            << std::setw(14) << "FBCX"
            << std::setw(14) << "FBCY" << std::endl;
    }

    for (const std::unique_ptr<Mesh>& msh : _meshs) {
        for (const Node& node : msh->nodes) {
            if (node.mn > 0.0) {
                nfile << std::setw(8) << node.nid
                    << std::setw(14) << std::scientific << node.mn
                    << std::setw(14) << node.xn[0]
                    << std::setw(14) << node.xn[1]
                    << std::setw(14) << node.vn[0]
                    << std::setw(14) << node.vn[1]
                    << std::setw(14) << node.an[0]
                    << std::setw(14) << node.an[1]
                    << std::setw(14) << node.fint[0]
                    << std::setw(14) << node.fint[1]
                    << std::setw(14) << node.fext[0]
                    << std::setw(14) << node.fext[1]
                    << std::setw(14) << node.fbc[0]
                    << std::setw(14) << node.fbc[1] << std::endl;
            }
        }
    }
    nfile.close();
    // std::cout << nfile_name << " has been output" << std::endl;
}

// ****************************    BOUNDARY    ***************************************
// Constructor
Boundary::Boundary(const std::string& form, const Vector2D& p1, const Vector2D& p2)
    : _form{ form }, _mu{ 0.0 }, _p1{ p1 }, _p2{ p2 }, _nbc{} 
{
    Vector2D diff = _p2 - _p1;

    // Assume _basey is the vertical vector of _basex
    _exl = diff.normalized();
    _eyl = Vector2D(-_exl[1], _exl[0]);
    //std::cout << _basex << ", " << _basey << std::endl;

    // calculate transform matrix (global -> local)
    double exlx = _exl.dot(ex);
    double exly = _exl.dot(ey);
    double eylx = _eyl.dot(ex);
    double eyly = _exl.dot(ey);

    Vector2D exlxy{ exlx, exly };
    Vector2D eylxy{ eylx, eyly };
    _Txlx.push_back(exlxy);
    _Txlx.push_back(eylxy);

    // calculate transform matrix (local -> global)
    Vector2D exylx{ exlx, eylx };
    Vector2D exyly{ exly, eyly };
    _Txxl.push_back(exylx);
    _Txxl.push_back(exyly);
}

// Destructor
Boundary::~Boundary() {}

// Setters
void Boundary::setForm(const std::string& form)
{
    _form = form;
}

void Boundary::setMu(double mu)
{
    _mu = mu;
}

void Boundary::setP1(const Vector2D& p1)
{
    _p1 = p1;
}

void Boundary::setP2(const Vector2D& p2)
{
    _p2 = p2;
}

void Boundary::setNbc(const Vector2D& nbc)
{
    _nbc = nbc;
}

// Getters (return by reference)
const std::string& Boundary::form() const
{
    return _form;
}

double Boundary::mu() const
{
    return _mu;
}

const Vector2D& Boundary::p1() const
{
    return _p1;
}

const Vector2D& Boundary::p2() const
{
    return _p2;
}

const Vector2D& Boundary::nbc() const
{
    return _nbc;
}

Vector2D Boundary::transformGlobal(const Vector2D& localPosition)
{
    Vector2D diff = localPosition - _p1;
    return Vector2D{ diff.dot(_Txxl[0]), diff.dot(_Txxl[1]) };
}

Vector2D Boundary::transformLocal(const Vector2D& globalPosition)
{
    Vector2D temp{ globalPosition.dot(_Txlx[0]), globalPosition.dot(_Txlx[1]) };
    return temp + _p1;
}

Vector2D Boundary::newtonMethod(const Vector2D& globalPosition)
{
    double xo = globalPosition[0];
    double yo = globalPosition[1];
    double x = xo;
    double y = 0.0;
    Vector2D Tt{ 1.0, 0.0 };

    int num{};
    while (true)
    {
        double A{};
        double dA{};
        double B{};
        double dB{};
        double N{};
        Vector2D n{};
        if (yo - y == 0.0)
        {
            n = _nbc;
            if (x < xo)
                n *= (-1.0);
        } 

        else
        {
            n = Vector2D(xo - x, yo - y) / std::sqrt((xo - x) * (xo - x) + (yo - y) * (yo - y));
            A = 1.0 / std::sqrt((xo - x) * (xo - x) + (yo - y) * (yo - y));
            dA = -0.5 * std::pow((xo - x) * (xo - x) + (yo - y) * (yo - y), -1.5) * (-2.0 * (xo - x) - 2.0 * (yo - y) * 0.0);
            B = Tt[0] * (xo - x) + Tt[1] * (yo - y);
            dB = -Tt[0] + 0.0 * (yo - y) + Tt[1] * 0.0;
        }

        N = Tt.dot(n);
        if (std::abs(N) <= TOL || num > 100)
            return Vector2D(x, y);

        double dN = dA * B + A * dB;
        if (dN != 0.0)
            x -= N / dN;
        y = 0.0;
        num++;
    }
}

// define layer function
static double layer(double r, const std::function<double(double)>& decayFunction)
{
    if (r >= UNITGRID)
        return 0.0;

    else if (0.0 < r && r < UNITGRID)
        return decayFunction(r);

    else
        return 1.0;
}
