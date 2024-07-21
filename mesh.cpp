#include "mesh.h"

// ****************************    MATERIAL    ***************************************
Material::Material(double rho, double K, double G) :
    rho{ rho }, K{ K }, G{ G }
{
    E = 9.0 * K * G / (3.0 * K + G);
    v = 0.5 * (3.0 * K - 2.0 * G) / (3.0 * K + G); // poisson

    // plane stress
    E1 = E / (1.0 - v * v);
    E2 = v * E / (1.0 - v * v);
}

Material::~Material() {}

void Material::verify_time_step(double dt)
{
    // P wave velocity->vel_p
    double M = K + 4.0 * G / 3.0;
    double vel_p = std::sqrt(M / rho);
    double dt_critical = UNITGRID / (10.0 * vel_p);
    if (dt > dt_critical)
    {
        std::cout << "dt = " << dt << std::endl;
        std::cout << "dt_critical = " << dt_critical << std::endl;
        std::cout << "need reset MPM time step" << std::endl;
        // Equivalent to Python's sys.exit()
        std::exit(EXIT_FAILURE);
    }
}

// ****************************    PARTICLE    ***************************************
// Constructor to initialize member variables
Particle::Particle() :
    pid{ 0 }, Vol{ 0.0 }, mp{ 0.0 }, xp{}, vp{}, Pp{}, bp{},
    ep{ 0.0, 0.0, 0.0 }, ssp{ 0.0, 0.0, 0.0 }, ineps{ 0.0, 0.0, 0.0 }, 
    N1{ 0.0 }, N2{ 0.0 }, N3{ 0.0 }, N4{ 0.0 },
    dN1{}, dN2{}, dN3{}, dN4{}
{
    // Initialization list is used to set initial values for member variables
}

// Destructor to clean up resources
Particle::~Particle() {}

void Particle::calculateMass(const Material& material)
{
    mp = material.rho * Vol;
}

void Particle::calculateMomentum()
{
    Pp = mp * vp;
}

void Particle::calculateSpecificStress(const Material& material)
{
    double sp[3]{};
    sp[0] = material.E1 * ep[0] + material.E2 * ep[1];
    sp[1] = material.E2 * ep[0] + material.E1 * ep[1];
    sp[2] = material.G * ep[2];

    ssp[0] = sp[0] / material.rho;
    ssp[1] = sp[1] / material.rho;
    ssp[2] = sp[2] / material.rho;
}

void Particle::show() const
{
    std::cout << "id : " << pid << std::endl;
    std::cout << "Volume : " << Vol << std::endl;
    std::cout << "mass : " << mp << std::endl;
    std::cout << "position : " << xp << std::endl;
    std::cout << "velocity : " << vp << std::endl;
    std::cout << "momentum : " << Pp << std::endl;
    std::cout << "specific stress : " << ssp[0] << ", " << ssp[1] << ", " << ssp[2] << std::endl;
}

// ****************************    NODE    ***************************************
// Constructor to initialize member variables
Node::Node() :
    nid{ 0 }, mn{ 0.0 }, xn{}, vn{}, an{}, pn{},
    fint{}, fext{}, ftot{}, fbc{}, bn{}, normal{}
{
    // Initialization list is used to set initial values for member variables
}

// Destructor to clean up resources
Node::~Node()
{
    // No dynamic memory to clean up for now
}

// ****************************    ELEMENT    ***************************************
Element::Element() :n1{ nullptr }, n2{ nullptr }, n3{ nullptr }, n4{ nullptr } {}

Element::~Element() {}

bool Element::contains(const Particle& particle) const
{
    // Use a simple bounding box check for demonstration purposes
    double minX = std::min({ n1->xn[0], n2->xn[0], n3->xn[0], n4->xn[0] });
    double maxX = std::max({ n1->xn[0], n2->xn[0], n3->xn[0], n4->xn[0] });
    double minY = std::min({ n1->xn[1], n2->xn[1], n3->xn[1], n4->xn[1] });
    double maxY = std::max({ n1->xn[1], n2->xn[1], n3->xn[1], n4->xn[1] });

    return particle.xp[0] >= minX && particle.xp[0] < maxX &&
           particle.xp[1] >= minY && particle.xp[1] < maxY;
}

// ****************************    MESH    ***************************************
Mesh::Mesh(const std::string& particleFile, const std::string& nodeFile, const Material& material)
    : material{ material }
{
    initParticleInfo(particleFile, material);
    initNodeInfo(nodeFile);
    createElements();
    createElementParticleMap();
}

Mesh::~Mesh() 
{
    // std::cout << "The mesh has been deleted\n";
}

void Mesh::initParticleInfo(const std::string& filePath, const Material& material)
{
    std::ifstream file(filePath);
    std::string line;
    while (std::getline(file, line))
    {
        // Skip lines that start with #
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        Particle ip{};
        iss >> ip.pid
            >> ip.Vol
            >> ip.xp[0] >> ip.xp[1]
            >> ip.vp[0] >> ip.vp[1]
            >> ip.ep[0] >> ip.ep[1] >> ip.ep[2]
            >> ip.bp[0] >> ip.bp[1];

        ip.calculateSpecificStress(material);

        if (ip.Vol != 0.0)
            particles.push_back(ip);
    }
    // std::cout << particles.size() << std::endl;
    file.close();
}

void Mesh::createElements()
{
    // Assume the grid is square
    int nodesPerSide = static_cast<int>(sqrt(nodes.size())) - 1;
    for (int i = 0; i < nodesPerSide; ++i)
    {
        for (int j = 0; j < nodesPerSide; ++j)
        {
            int n1_id = i * (nodesPerSide + 1) + j;
            int n2_id = n1_id + 1;
            int n3_id = (i + 1) * (nodesPerSide + 1) + j;
            int n4_id = n3_id + 1;

            Element e{};
            e.n1 = &nodes[n1_id];
            e.n2 = &nodes[n2_id];
            e.n3 = &nodes[n4_id];
            e.n4 = &nodes[n3_id];
            elements.push_back(std::move(e));
        }
    }
}

void Mesh::createElementParticleMap()
{
    for (int i = 0; i < particles.size(); ++i)
    {
        for (int j = 0; j < elements.size(); ++j)
        {
            if (elements[j].contains(particles[i]))
            {
                pem[i] = j;
                break;
            }
        }
    }
}

void Mesh::clearMapId()
{
    pem.clear();
}

//Element* Mesh::findElementForParticle(const Particle& particle)
// {
//    auto it = std::find_if(particles.begin(), particles.end(), [&](const Particle& p) { 
//        return p.xp[0] == particle.xp[0] && p.xp[1] == particle.xp[1]; 
//    });
// 
//    if (it != particles.end())
//    {
//        int particleIndex = std::distance(particles.begin(), it);
//        if (pem.find(particleIndex) != pem.end()) {
//            return &elements[pem[particleIndex]];
//        }
//    }
//    return nullptr;
//}

void Mesh::initNodeInfo(const std::string& filename)
{
    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    std::string line;
    int nodeCount;
    // Read the total number of nodes
    if (!(infile >> nodeCount))
    {
        std::cerr << "Error reading the number of nodes from file: " << filename << std::endl;
        return;
    }

    // Skip the next few lines
    for (int i = 0; i < 4; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Error reading header lines from file: " << filename << std::endl;
            return;
        }
    }

    // Read node data
    int id = 0;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        Node n{};
        if (iss >> n.xn[0] >> n.xn[1]) {
            n.nid = id++;
            nodes.emplace_back(n);
        }
    }
    infile.close();
}

void Mesh::showParticleInitInfo() const
{
    for (const Particle& p : particles)
    {
        std::cout << "Particle ID: " << p.pid << "\n"
                  << "Volume: " << p.Vol << "\n"
                  << "Position: (" << p.xp[0] << ", " << p.xp[1] << ")\n"
                  << "Velocity: (" << p.vp[0] << ", " << p.vp[1] << ")\n"
                  << "Momentum: (" << p.Pp[0] << ", " << p.Pp[1] << ")\n"
                  << "Body Force: (" << p.bp[0] << ", " << p.bp[1] << ")\n"
                  << "Strain: (" << p.ep[0] << ", " << p.ep[1] << ", " << p.ep[2] << ")\n"
                  << std::endl;
    }
}

void Mesh::showNodeInitInfo() const
{
    for (const Node& node : nodes)
    {
        std::cout << "Node ID: " << node.nid 
                  << ", Position: (" << node.xn[0] << ", " << node.xn[1] << ")\n";
    }
}

void Mesh::showElementInfo() const
{
    int eid = 0;
    for (const Element& elem : elements)
    {
        std::cout << "Element " << eid++
                  << " -node: "
                  << elem.n1->nid << " (" << elem.n1->xn[0] << ", " << elem.n1->xn[1] << "), "
                  << elem.n2->nid << " (" << elem.n2->xn[0] << ", " << elem.n2->xn[1] << "), "
                  << elem.n3->nid << " (" << elem.n3->xn[0] << ", " << elem.n3->xn[1] << "), "
                  << elem.n4->nid << " (" << elem.n4->xn[0] << ", " << elem.n4->xn[1] << ")\n";
    }
}

static double stable(double fij)
{
    if (fij < N_min)
        return N_min;

    else if (fij > (1.0 - N_min))
        return 1.0 - N_min;

    else
        return fij;
}


static double  stabledN(double dfij, double fij)
{
    double Le = 0.5 - std::abs(fij - 0.5);
    if (Le < L_cut)
        return fij * dfij * Le / L_cut;

    return fij * dfij;
}

double shapeN1(const double xi, const double eta) {
    return stable(0.5 * (1.0 - xi)) * stable(0.5 * (1.0 - eta));
}
double shapeN2(const double xi, const double eta) {
    return stable(0.5 * (1.0 + xi)) * stable(0.5 * (1.0 - eta));
}
double shapeN3(const double xi, const double eta) {
    return stable(0.5 * (1.0 + xi)) * stable(0.5 * (1.0 + eta));
}
double shapeN4(const double xi, const double eta) {
    return stable(0.5 * (1.0 - xi)) * stable(0.5 * (1.0 + eta));
}
double shapedN1x(const double eta) {
    return stabledN((0.5 * (-1.0) * dxi), stable(0.5 * (1.0 - eta)));
}
double shapedN1y(const double xi) {
    return stabledN((0.5 * (-1.0) * deta), stable(0.5 * (1.0 - xi)));
}
double shapedN2x(const double eta) {
    return stabledN((0.5 * (1.0) * dxi), stable(0.5 * (1.0 - eta)));
}
double shapedN2y(const double xi) {
    return stabledN((0.5 * (-1.0) * deta), stable(0.5 * (1.0 + xi)));
}
double shapedN3x(const double eta) {
    return stabledN((0.5 * (1.0) * dxi), stable(0.5 * (1.0 + eta)));
}
double shapedN3y(const double xi) {
    return stabledN((0.5 * (1.0) * deta), stable(0.5 * (1.0 + xi)));
}
double shapedN4x(const double eta) {
    return stabledN((0.5 * (-1.0) * dxi), stable(0.5 * (1.0 + eta)));
}
double shapedN4y(const double xi) {
    return stabledN((0.5 * (1.0) * deta), stable(0.5 * (1.0 - xi)));
}
    