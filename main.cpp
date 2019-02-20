#include <vector>
#include <algorithm>
#include <memory>
#include <iterator>
#include <string>
#include <iostream>
#include <fstream>

int NCONTAINERS = 4500;

class Clients
{
public:
    std::vector<double> M_freq;
    std::vector<int> M_locations;

    Clients(std::string const& filename );
};

Clients::Clients(std::string const& filename )
{
    std::ifstream myfile;
    myfile.open(filename);
    if( myfile.is_open() )
    {
        std::string line;
        while( getline(myfile,line) )
        {
            int oldPos = 0;
            std::string delim = ",";
            auto pos = line.find(delim);
            while( pos != std::string::npos )
            {
                std::string part = line.substr(oldPos,pos-oldPos);
                auto pospc = part.find('%');
                if( pospc != std::string::npos )
                {
                    M_freq.push_back(std::stod(part.substr(1,pospc)));
                }
                oldPos = pos;
                pos = line.find(delim,pos+1);
            }
        }
        myfile.close();
        double fq = 0;
        std::sort(M_freq.begin(), M_freq.end());
        std::reverse(M_freq.begin(), M_freq.end());
        for( auto const& f : M_freq )
        {
            M_locations.push_back(int(NCONTAINERS/5*f/100.));
            fq += f;
        }
    }
    else
        std::cout << "unable to open file " << filename << std::endl;
}

class Bay;
class Block
{
    using block_ptrtype = std::shared_ptr<Block>;

public:
    std::string M_zone;
    int M_number;
    int M_locations;
    double M_weight;
    std::shared_ptr<Bay> M_in;
    std::shared_ptr<Bay> M_out;

public:
    Block(std::string const& zone, int number, int locations, std::shared_ptr<Bay> const& in, std::shared_ptr<Bay> const& out, double weight = 0);
    void setWeight(double weight) { M_weight = weight; }
    int number() const { return M_number; }
    int locations() const { return M_locations; }
    int locationsWeighted() const { return M_weight == 0 ? 0 : M_locations; }
    double weight() const { return M_weight; }
    double distance() const;
    double cost() const { return distance()*M_weight*M_locations; }
};

class Bay
{
    using block_ptrtype = std::shared_ptr<Block>;
public:
    double M_distanceIn;
    double M_distanceOut;
    std::vector<block_ptrtype> M_blocksIn;
    std::vector<block_ptrtype> M_blocksOut;

public:
    Bay(double distanceIn, double distanceOut, std::vector<block_ptrtype> in, std::vector<block_ptrtype> out);
    void setBlockIn(std::vector<block_ptrtype> in) { M_blocksIn = in; }
    void setBlockOut(std::vector<block_ptrtype> out) { M_blocksOut = out; }
    void addBlockIn(block_ptrtype in) { M_blocksIn.push_back(in); }
    void addBlockOut(block_ptrtype out) { M_blocksOut.push_back(out); }
    void addBlockInFront(block_ptrtype in) { M_blocksIn.insert(M_blocksIn.begin(),in); }
    void addBlockOutFront(block_ptrtype out) { M_blocksOut.insert(M_blocksOut.begin(),out); }
    double distanceIn() const { return M_distanceIn; }
    double distanceOut() const { return M_distanceOut; }
    int lengthIn() const { return M_blocksIn.size()*20; }
    int lengthOut() const { return M_blocksOut.size()*20; }
    // int fromIn() const { return M_blocksIn.front()->number(); }
    // int fromOut() const { return M_blocksOut.front()->number(); }
    // int toIn() const { return M_blocksIn.back()->number(); }
    // int toOut() const { return M_blocksOut.back()->number(); }
    int nbLocationsIn() const;
    int nbLocationsWeightedIn() const;
};

class Port
{
    using block_ptrtype = std::shared_ptr<Block>;
    using bay_ptrtype = std::shared_ptr<Bay>;
    
    std::vector<bay_ptrtype> M_bays;
    std::vector<block_ptrtype> M_blocks;

public:
    Port() {};
    static Port generateCurrentPAS();
    static Port generateRomain();
    void setBays(std::vector<bay_ptrtype> const& bays) { M_bays = bays; }
    void setBlocks(std::vector<block_ptrtype> const& blocks) { M_blocks = blocks; }
    double cost() const;
    int nbLocations() const;
    int nbLocationsWeighted() const;
    int nbContainers() const { return 5*nbLocations(); }
    int nbContainersWeighted() const { return 5*nbLocationsWeighted(); }
    int nbBlocks() const { return M_blocks.size(); }
    int nbBlocksWeighted() const;
    void setWeight( std::vector<int> nbLocations, std::vector<double> weigths);
};

Block::Block(std::string const& zone, int number, int locations, std::shared_ptr<Bay> const& in, std::shared_ptr<Bay> const& out, double weight)
    : M_zone(zone),
      M_number(number),
      M_locations(locations),
      M_weight(weight),
      M_in(in),
      M_out(out)
{}

double Block::distance() const
{
    double d = M_in->distanceIn() + M_out->distanceOut();    
    auto itIn = std::find_if(M_in->M_blocksIn.begin(), M_in->M_blocksIn.end(), [&](block_ptrtype const& b) { return b->M_zone == this->M_zone && b->M_number == this->M_number; } );;
    d += std::distance(M_in->M_blocksIn.begin(), itIn);
    auto itOut = std::find_if(M_out->M_blocksOut.begin(), M_out->M_blocksOut.end(), [&](block_ptrtype const& b) { return b->M_zone == this->M_zone && b->M_number == this->M_number; } );;
    d += std::distance(M_out->M_blocksOut.begin(), itOut);
    return d;
}

Bay::Bay(double distanceIn, double distanceOut, std::vector<block_ptrtype> in = std::vector<block_ptrtype>(), std::vector<block_ptrtype> out = std::vector<block_ptrtype>())
    : M_distanceIn(distanceIn),
      M_distanceOut(distanceOut),
      M_blocksIn(in),
      M_blocksOut(out)
{}

int Bay::nbLocationsIn() const
{
    int nb = 0;
    for(auto const& bl : M_blocksIn )
        nb += bl->locations();
    return nb;
}

int Bay::nbLocationsWeightedIn() const
{
    int nb = 0;
    for(auto const& bl : M_blocksIn )
        nb += bl->locationsWeighted();
    return nb;
}

Port Port::generateCurrentPAS()
{
    auto p = Port();
    std::vector<bay_ptrtype> bays;
    bays.push_back(std::make_shared<Bay>(464,208));
    bays.push_back(std::make_shared<Bay>(352,96));
    bays.push_back(std::make_shared<Bay>(256,0));
    bays.push_back(std::make_shared<Bay>(160,96));
    bays.push_back(std::make_shared<Bay>(0,256));
    bays.push_back(std::make_shared<Bay>(80,336));
    bays.push_back(std::make_shared<Bay>(160,416));
    bays.push_back(std::make_shared<Bay>(320,576));
    bays.push_back(std::make_shared<Bay>(400,656));
    bays.push_back(std::make_shared<Bay>(480,736));
    bays.push_back(std::make_shared<Bay>(560,816));
    bays.push_back(std::make_shared<Bay>(640,896));
    bays.push_back(std::make_shared<Bay>(720,976));

    std::vector<block_ptrtype> blocks;
    for(int i = 3; i < 13; ++i)
    {
        blocks.push_back(std::make_shared<Block>("H2",i,10,bays[1],bays[0]));
        bays[0]->addBlockOut(blocks.back());
        bays[1]->addBlockInFront(blocks.back());
    }
    for(int i = 1; i < 15; ++i)
    {
        blocks.push_back(std::make_shared<Block>("H5",i,6,bays[2],bays[1]));
        bays[2]->addBlockIn(blocks.back());
        bays[1]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("H7",i,6,bays[3],bays[2]));
        bays[3]->addBlockIn(blocks.back());
        bays[2]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("I1",i,14,bays[4],bays[3]));
        bays[4]->addBlockIn(blocks.back());
        bays[3]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("I5",i,4,bays[5],bays[4]));
        bays[5]->addBlockIn(blocks.back());
        bays[4]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("I7",i,4,bays[6],bays[5]));
        bays[6]->addBlockIn(blocks.back());
        bays[5]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("J1",i,14,bays[7],bays[6]));
        bays[7]->addBlockIn(blocks.back());
        bays[6]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("J5",i,4,bays[8],bays[7]));
        bays[8]->addBlockIn(blocks.back());
        bays[7]->addBlockOut(blocks.back());
        blocks.push_back(std::make_shared<Block>("J7",i,4,bays[9],bays[8]));
        bays[9]->addBlockIn(blocks.back());
        bays[8]->addBlockOut(blocks.back());
    }
    for(int i = 5; i < 15; ++i)
    {
        blocks.push_back(std::make_shared<Block>("K1",i,4,bays[10],bays[9]));
        bays[10]->addBlockIn(blocks.back());
        bays[9]->addBlockOut(blocks.back());
    }
    for(int i = 8; i < 15; ++i)
    {
        blocks.push_back(std::make_shared<Block>("K3",i,4,bays[11],bays[10]));
        bays[11]->addBlockIn(blocks.back());
        bays[10]->addBlockOut(blocks.back());
    }
    blocks.push_back(std::make_shared<Block>("K5",10,2,bays[12],bays[11]));
    bays[12]->addBlockIn(blocks.back());
    bays[11]->addBlockOut(blocks.back());
    for(int i = 11; i < 15; ++i)
    {
        blocks.push_back(std::make_shared<Block>("K5",i,4,bays[12],bays[11]));
        bays[12]->addBlockIn(blocks.back());
        bays[11]->addBlockOut(blocks.back());
    }

    p.setBays(bays);
    p.setBlocks(blocks);
    return p;
}

Port Port::generateRomain()
{
    auto p = Port();
    std::vector<bay_ptrtype> bays;
    bays.push_back(std::make_shared<Bay>(96,328));
    bays.push_back(std::make_shared<Bay>(0,328));
    bays.push_back(std::make_shared<Bay>(452,876));
    bays.push_back(std::make_shared<Bay>(456,568));
    bays.push_back(std::make_shared<Bay>(540,788));
    bays.push_back(std::make_shared<Bay>(540,788));
    bays.push_back(std::make_shared<Bay>(264,160));
    bays.push_back(std::make_shared<Bay>(264,160));
    bays.push_back(std::make_shared<Bay>(620,708));
    bays.push_back(std::make_shared<Bay>(712,448));
    bays.push_back(std::make_shared<Bay>(700,628));
    bays.push_back(std::make_shared<Bay>(424,0));
    bays.push_back(std::make_shared<Bay>(424,0));

    std::vector<block_ptrtype> blocks;
    for(int i = 0; i < 10; ++i)
    {
        blocks.push_back(std::make_shared<Block>("A1",i,15,bays[0],bays[6]));
        bays[0]->addBlockIn(blocks.back());
        bays[6]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 6; ++i)
    {
        blocks.push_back(std::make_shared<Block>("A2",i,5,bays[0],bays[3]));
        bays[0]->addBlockIn(blocks.back());
        bays[3]->addBlockOutFront(blocks.back());
    }
    for(int i = 0; i < 6; ++i)
    {
        blocks.push_back(std::make_shared<Block>("A3",i,4,bays[3],bays[6]));
        bays[3]->addBlockInFront(blocks.back());
        bays[6]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 15; ++i)
    {
        blocks.push_back(std::make_shared<Block>("B1",i,15,bays[1],bays[7]));
        bays[1]->addBlockIn(blocks.back());
        bays[7]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 10; ++i)
    {
        blocks.push_back(std::make_shared<Block>("B2",i,5,bays[1],bays[4]));
        bays[1]->addBlockIn(blocks.back());
        bays[4]->addBlockOutFront(blocks.back());
    }
    for(int i = 0; i < 10; ++i)
    {
        blocks.push_back(std::make_shared<Block>("B3",i,4,bays[4],bays[7]));
        bays[4]->addBlockInFront(blocks.back());
        bays[7]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 25; ++i)
    {
        blocks.push_back(std::make_shared<Block>("C1",i,5,bays[2],bays[5]));
        bays[2]->addBlockIn(blocks.back());
        bays[5]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 25; ++i)
    {
        blocks.push_back(std::make_shared<Block>("C2",i,4,bays[5],bays[8]));
        bays[5]->addBlockIn(blocks.back());
        bays[8]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 10; ++i)
    {
        blocks.push_back(std::make_shared<Block>("D1",i,14,bays[6],bays[11]));
        bays[6]->addBlockIn(blocks.back());
        bays[11]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 6; ++i)
    {
        blocks.push_back(std::make_shared<Block>("D2",i,4,bays[6],bays[9]));
        bays[6]->addBlockIn(blocks.back());
        bays[9]->addBlockOutFront(blocks.back());
    }
    for(int i = 0; i < 6; ++i)
    {
        blocks.push_back(std::make_shared<Block>("D3",i,4,bays[9],bays[11]));
        bays[9]->addBlockInFront(blocks.back());
        bays[11]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 15; ++i)
    {
        blocks.push_back(std::make_shared<Block>("E1",i,14,bays[7],bays[12]));
        bays[7]->addBlockIn(blocks.back());
        bays[12]->addBlockOut(blocks.back());
    }
    for(int i = 0; i < 10; ++i)
    {
        blocks.push_back(std::make_shared<Block>("E2",i,4,bays[7],bays[10]));
        bays[7]->addBlockIn(blocks.back());
        bays[10]->addBlockOutFront(blocks.back());
    }
    for(int i = 0; i < 10; ++i)
    {
        blocks.push_back(std::make_shared<Block>("E3",i,4,bays[10],bays[12]));
        bays[10]->addBlockInFront(blocks.back());
        bays[12]->addBlockOut(blocks.back());
    }

    p.setBays(bays);
    p.setBlocks(blocks);
    return p;
}

int Port::nbLocations() const
{
    int nb = 0;
    for( auto const& ba : M_bays )
        nb += ba->nbLocationsIn();
    return nb;
}

int Port::nbLocationsWeighted() const
{
    int nb = 0;
    for( auto const& ba : M_bays )
        nb += ba->nbLocationsWeightedIn();
    return nb;
}

int Port::nbBlocksWeighted() const
{
    int nb = 0;
    for( auto const& bl : M_blocks )
        nb += bl->weight() != 0 ? 1 : 0;
    return nb;
}

void Port::setWeight( std::vector<int> nbLocations, std::vector<double> weights)
{
    std::sort(M_blocks.begin(), M_blocks.end(), [&](block_ptrtype const& b1, block_ptrtype const& b2) {return b1->distance() < b2->distance(); } );
    auto it = M_blocks.begin();
    for( int i = 0; i < nbLocations.size() && i < 50; ++i)
    {
        int Xc = nbLocations[i];
        while( nbLocations[i] > 0 && it != M_blocks.end() )
        {
            (*it)->setWeight(weights[i]/Xc);
            nbLocations[i] -= (*it)->locations();
            std::advance(it,1);
        }
    }
}

double Port::cost() const
{
    double c = 0;
    for( auto t : M_blocks)
        c += t->cost();

    return c/this->nbBlocks();
}

int main(int argc, char** argv)
{
    auto clients = Clients("clientlist.csv");
    Port p;
    if( argc == 1 || atoi(argv[1]) == 1)
        p = Port::generateCurrentPAS();
    else if( atoi(argv[1]) == 2 )
        p = Port::generateRomain();
    else
        p = Port::generateCurrentPAS();
    // auto p = Port::generateCurrentPAS();
    // p.setWeight({{100,80,60,40,34,34,34,34,34,34,12,12,12,12,12,12,12,12,12,12}},{{0.15,0.12,0.11,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.022}});
    p.setWeight(clients.M_locations, clients.M_freq);
    auto c = p.cost();
    std::cout << "port:"
              << "\n- conteneurs: capacite totale = " << p.nbContainers() << " | occupe = " << p.nbContainersWeighted()
              << "\n- emplacements: capacite totale = " << p.nbLocations() << " | occupe = " << p.nbLocationsWeighted()
              << "\n- travees: capacite totale = " << p.nbBlocks() << " | occupe = " << p.nbBlocksWeighted()
              << "\n\ncout = " << c << std::endl;

    return 1;
}
