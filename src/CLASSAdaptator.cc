#include "CLASSAdaptator.h"


using cyclus::Nuc;
using cyclus::Material;
using cyclus::Composition;

//#define cyDBGL		std::cout << __FILE__ << " : " << __LINE__ << " [" << __FUNCTION__ << "]" << std::endl;
#define cyDBGL ;

namespace cyclass {


    CLASSAdaptator::CLASSAdaptator(std::string EqModel, std::string XSModel, std::string BatemanSolver){

        cyDBGL
        myPhysicsModel = new PhysicsModels();
        if( EqModel == "xx"){
            myPhysicsModel->SetEquivlalenceModel(new EquivalenceModel());
        }
        else if(EqModel == "xx" ){
            myPhysicsModel->SetEquivlalenceModel(new EquivalenceModel());
        }

        if( XSModel == "xx"){
            myPhysicsModel->SetXSModel(new XSModel());
        }
        else if(EqModel == "xx" ){
            myPhysicsModel->SetXSModel(new XSModel());
        }


        if( BatemanSolver == "ii"){
            myPhysicsModel->SetIrradiationModel(new IrradiationModel());
        }
        else if(BatemanSolver == "ii" ){
            myPhysicsModel->SetIrradiationModel(new IrradiationModel();
        }

        cyDBGL

    }

    //________________________________________________________________________
    float CLASSAdaptator::GetEnrichment(cyclus::Composition::Ptr c_fissil,
                                        cyclus::Composition::Ptr c_fertil,
                                        double BurnUp){
        cyDBGL

        double val = 0;

        IsotopicVector IV_fissil = CYCLUS2CLASS(c_fissil);
        IsotopicVector IV_fertil = CYCLUS2CLASS(c_fertil);

        val= myPhysicsModel->GetEquivalenceModel()->GetFissileMolarFraction(IV_fissil, IV_fertil, BurnUp);

        cyDBGL

        return val; //
    }


    //________________________________________________________________________
    float MLPBUsolver::GetBU(cyclus::Composition::Ptr fuel, double eps)    {
        cyDBGL

        float BU_max = 80;
        float BU_min = 10;
        cyDBGL

        cyclus::Composition::Ptr fuel_fissil = ExtractAccordinglist( fuel, fissil_list);
        cyclus::Composition::Ptr fuel_fertil = ExtractAccordinglist( fuel, fertil_list);



        if( std::abs(AtomIn(fuel_fertil) + AtomIn(fuel_fissil) - AtomIn(fuel)) > 1e-10 ){
            std::cout << "You fuel has nuclei that this model could not manage.."<< std::endl;
            exit(1);
        }


        float rho_target = AtomIn(fuel_fissil)/AtomIn(fuel);
        float rho_min = GetEnrichment(fuel_fissil, fuel_fertil, BU_min);
        float rho_max = GetEnrichment(fuel_fissil, fuel_fertil, BU_max);
        float BU_estimation = BU_max;
        float rho_estimated = rho_max;

        do {
            //Update BU_estimation

            if (rho_estimated > rho_target){
                rho_max = rho_estimated;
                BU_max = BU_estimation;

            } else {
                rho_min = rho_estimated;
                BU_min = BU_estimation;
            }

            BU_estimation = (BU_max+BU_min)/2.;

            rho_estimated = GetEnrichment(fuel_fissil, fuel_fertil, BU_estimation);


        }while( std::abs(rho_target - rho_estimated)/(rho_target/2.+rho_estimated/2.) > eps );

        cyDBGL
        return BU_estimation;

    }


    //________________________________________________________________________
    //________________________________________________________________________
    //________________________________________________________________________
    double AtomIn(CompMap Source){

        double total = 0;

        CompMap::iterator it;

        for(it = Source.begin(); it != Source.end(); it++)
            total += it->second;

        return total;
    }

    //________________________________________________________________________
    double AtomIn(cyclus::Composition::Ptr Source) { return AtomIn(Source->atom());};

    //________________________________________________________________________
    cyclus::Composition::Ptr ExtractAccordinglist(
                                                  cyclus::Composition::Ptr source,
                                                  cyclus::Composition::Ptr list){

        //create Output Composition
        CompMap separatedCompo;

        // Extract Nuc map from source compo & list...
        CompMap sourceComp = source->atom();
        CompMap ListComp = list->atom();

        CompMap::iterator it;

        // Fill output composition
        for (it = ListComp.begin(); it != ListComp.end(); it++) {
            CompMap::iterator it2 = sourceComp.find( it->first );

            if(it2 != sourceComp.end())
                separatedCompo.insert( std::pair<Nuc,double>(it->first, it2->second) );

        }

        return Composition::CreateFromAtom(separatedCompo);
    }

    //________________________________________________________________________
    CompMap NormalizeComp( CompMap source, double norm ){
        double Total = AtomIn(source);
        CompMap::iterator it;
        for(it = source.begin(); it != source.end(); it++)
            it->second *= norm/Total;

        return source;
    }

    //________________________________________________________________________
    CompMap operator+( CompMap const& IVa, CompMap const& IVb){

        CompMap IVtmp = IVa;
        CompMap IVbtmp = IVb;
        CompMap::iterator it;

        for(it = IVbtmp.begin(); it != IVbtmp.end(); it++){

            std::pair<CompMap::iterator, bool> IResult;
            IResult = IVtmp.insert(std::pair<Nuc, double> (it->first, it->second));
            if(!IResult.second)
                IResult.first->second += it->second;

        }
        return IVtmp;
    }


    //________________________________________________________________________
    CompMap operator*(CompMap const& IVA, double F){

        CompMap IVtmp = IVA;
        CompMap::iterator it;

        for(it = IVtmp.begin(); it != IVtmp.end(); it++)
            it->second *= F;

        return IVtmp;
    }

    //________________________________________________________________________
    void Print(cyclus::Composition::Ptr compo ){
        Print(compo->atom());
    }

    //________________________________________________________________________
    void Print(CompMap compo){
        CompMap::iterator it;
        
        for (it = compo.begin(); it != compo.end(); it++){
            std::cout << it->first << " " << it->second << std::endl;
        }
    }
    
    //________________________________________________________________________
    CompMap operator-(CompMap const& IVa, CompMap const& IVb) { return IVa + (-1*IVb); };

    //________________________________________________________________________
    CompMap operator/(CompMap const& IVA, double F) { return IVA * (1/F); };

    //________________________________________________________________________
    CompMap operator*(double F, CompMap const& IVA) { return IVA*F; };

    //________________________________________________________________________
    IsotopicVector CYCLUS2CLASS(CompMap c_compo){

        IsotopicVector IV;
        CompMap::iterator it;

        for( it = c_compo.begin(); it != c_compo.end(); it++){
            int Z = (int)c_fissil->first /10000000 ;
            int A = (int)(c_fissil->first - Z)/10000;
            int I = 0;
            if( c_fissil->first - Z - A > 0 ) I++;
            IV += c_fissil->second() * ZAI(Z, A, I);
        }

        return IV;
    }

    //________________________________________________________________________
    CompMap CLASS2CYCLUS(IsotopicVector IV){

        CompoMap myCompoMap;

        map<ZAI, double> IsotopicMap = IV.GetIsotopicQuantity();
        map<ZAI, double>::iterator it;

        for(it = IsotopicMap.begin(); it != IsotopicMap.end(); it++){

            Nuc MyNuc = it->first.Z()*10000000 + it->first.A()*10000;
            if( I >0){
                std::cout << "you are screwed!! Isomer conversion from CLASS to CYCLUS... did not work yet..." <<std::endl;
                exit(1);
            }

            myCompMap.insert( pair<Muc, double> (MyNuc, it->second) );
        }
        return myCompoMap;


    }

    
    //________________________________________________________________________
    
    
    
    
    
} // namespace cyclass
