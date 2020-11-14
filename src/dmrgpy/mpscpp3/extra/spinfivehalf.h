//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_SPINFIVEHALF_H
#define __ITENSOR_SPINFIVEHALF_H
#include "itensor/mps/siteset.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/str.h"

namespace itensor {

class SpinFiveHalf : public SiteSet
    {
    public:

    SpinFiveHalf() { }

    SpinFiveHalf(int N, 
            Args const& args = Args::global());

    SpinFiveHalf(std::vector<Index> const& inds);

    void
    read(std::istream& s);

    };


class SpinFiveHalfSite
    {
    Index s;
    public:

    SpinFiveHalfSite(Index I) : s(I) { }

    SpinFiveHalfSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,S=3/2");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        auto conserveqns = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveqns);
        if(conserveSz)
            {
            s = Index(QN({"Sz",+5}),1,
                      QN({"Sz",+3}),1,
                      QN({"Sz",+1}),1,
                      QN({"Sz",-1}),1,
                      QN({"Sz",-3}),1,
                      QN({"Sz",-5}),1,Out,ts);
            }
        else
            {
            s = Index(4,ts);
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
	if (state == "Up" || state == "5")
            {
            return s(1);
            }
        else if (state == "Upi" || state == "3")
            {
            return s(2);
            }
        else if (state == "Upii" || state == "1")
            {
            return s(3);
            }
        else if (state == "Dnii" || state == "-1")
            {
            return s(4);
            }
        else if (state == "Dni" || state == "-3")
            {
            return s(5);
            }
        else if (state == "Dn" || state == "-5")
            {
            return s(6);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {

        const Real val1 = std::sqrt(5.0)/2.0;
        const Real val2 = std::sqrt(8.0)/2.0;
        const Real val3 = std::sqrt(9.0)/2.0;

        auto sP = prime(s);

        auto Up  = s(1);
        auto UpP = sP(1);
        auto Upi = s(2);
        auto UpiP = sP(2);
        auto Upii = s(3);
        auto UpiiP = sP(3);
        auto Dnii = s(4);
        auto DniiP = sP(4);
        auto Dni = s(5);
        auto DniP = sP(5);
        auto Dn  = s(6);
        auto DnP = sP(6);


        auto Op = ITensor(dag(s),sP);

        if(opname == "Sz")
            {
            Op.set(Up,UpP,+2.5);
            Op.set(Upi,UpiP,+1.5);
            Op.set(Upii,UpiiP,+0.5);
            Op.set(Dnii,DniiP,-0.5);
            Op.set(Dni,DniP,-1.5);
            Op.set(Dn,DnP,-2.5);
            }
        else
        if(opname == "Sx")
            {
            Op.set(Up,UpiP,val1);
            Op.set(Upi,UpP,val1);
            Op.set(Upi,UpiiP,val2);
            Op.set(Upii,UpiP,val2);
            Op.set(Upii,DniiP,val3);
            Op.set(Dnii,UpiiP,val3);
            Op.set(Dnii,DniP,val2);
            Op.set(Dni,DniiP,val2);
            Op.set(Dni,DnP,val1);
            Op.set(Dn,DniP,val1);
            }
        else
        if(opname == "Sy")
            {
            Op.set(Up,UpiP,val1*Cplx_i);
            Op.set(Upi,UpP,-val1*Cplx_i);
            Op.set(Upi,UpiiP,val2*Cplx_i);
            Op.set(Upii,UpiP,-val2*Cplx_i);
            Op.set(Upii,DniiP,val3*Cplx_i);
            Op.set(Dnii,UpiiP,-val3*Cplx_i);
            Op.set(Dnii,DniP,val2*Cplx_i);
            Op.set(Dni,DniiP,-val2*Cplx_i);
            Op.set(Dni,DnP,val1*Cplx_i);
            Op.set(Dn,DniP,-val1*Cplx_i);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    SpinFiveHalfSite(int n, Args const& args = Args::global())
        {
        *this = SpinFiveHalfSite({args,"SiteNumber=",n});
        }

    };

inline SpinFiveHalf::
SpinFiveHalf(std::vector<Index> const& inds)
    {
    int N = inds.size();
    auto sites = SiteStore(N);
    for(int j = 1, i = 0; j <= N; ++j, ++i)
        {
        auto& Ii = inds.at(i);
        if(dim(Ii) != 3)
            {
            printfln("Index at entry %d = %s",i,Ii);
            Error("Only S=1 IQIndices allowed in SpinFiveHalf(vector<Index>) constructor");
            }
        sites.set(j,SpinFiveHalfSite(Ii));
        }
    SiteSet::init(std::move(sites));
    }

inline SpinFiveHalf::
SpinFiveHalf(int N, 
        Args const& args)
    {
    auto shedge = args.getBool("SHalfEdge",false);
    auto Lshedge = args.getBool("SHalfLeftEdge",false);

    auto sites = SiteStore(N);

    auto start = 1;
    if(shedge || Lshedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site 1");
        sites.set(1,SpinHalfSite(1,args));
        start = 2;
        }

    for(int j = start; j < N; ++j)
        {
        sites.set(j,SpinFiveHalfSite(j,args));
        }

    if(shedge)
        {
        if(args.getBool("Verbose",false)) println("Placing a S=1/2 at site N=",N);
        sites.set(N,SpinHalfSite(N,args));
        }
    else
        {
        sites.set(N,SpinFiveHalfSite(N,args));
        }

    SiteSet::init(std::move(sites));
    }

void inline SpinFiveHalf::
read(std::istream& s)
    {
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        for(int j = 1; j <= N; ++j) 
            {
            auto I = Index{};
            I.read(s);
            if(dim(I) == 3) store.set(j,SpinFiveHalfSite(I));
            else if(dim(I) == 2) store.set(j,SpinHalfSite(I));
            else Error(format("SpinFiveHalf cannot read index of size %d",dim(I)));
            }
        init(std::move(store));
        }
    }

} //namespace itensor

#endif
