
#include "denseMBB.h"

//denseMBB::denseMBB( bgraph& g) {

//    bg = g;
//    maxsize = 0;


//}



//in general Ca shall be smaller than deg
void denseMBB::startEnu() {

    unordered_set<int> A;
    unordered_set<int> Ca;

    unordered_set<int> B;
    unordered_set<int> Cb;

    //L
    int v = 0;
    A.emplace(v);

    for(int i=1;i<bg.adjL.size();i++){
        Ca.emplace(i);

    }
    unordered_set<int> newCb;

    //R
    for(int i=1;i<bg.adjR.size();i++){
        Cb.emplace(i);
    }
    vector<int> nov = bg.adjL[v];
    for (auto &v:nov){
        newCb.emplace(v);
    }


    cout<<"alg starts"<<endl;

    recEnu(B,newCb,A,Ca,1);


    A.erase(v);
    recEnu(A,Ca,B,Cb,0);
}


void denseMBB::recEnu(unordered_set<int> &A, unordered_set<int> &Ca,
                      unordered_set<int> &B, unordered_set<int> &Cb,
                      int  flag) {

    //check pruning conditions
    //cou++;
    //cout<<""<<cou<<endl;
    if(min((A.size()+Ca.size()),(B.size() + Cb.size()))<=maxsize){
        //cout<<"pruning 1"<<endl;
        return;
    }
    //cout<<"test time 1"<<endl;
    unordered_set<int> newA, newB, newCa, newCb; //using memory in stack, which is efficient
    newA.insert(A.begin(),A.end());
    newB.insert(B.begin(),B.end());
    newCa.insert(Ca.begin(),Ca.end());
    newCb.insert(Cb.begin(),Cb.end());
    //cout<<"test time 2"<<endl;
    if(newCa.size() == 0 || newCb.size() == 0) {
        int currentsize = newA.size() + newCa.size();
        if (currentsize > newB.size() + newCb.size()) {
            currentsize = newB.size() + newCb.size();
        }

        //cout<<end<<endl;
        if (currentsize > maxsize) {
            maxsize = currentsize;

            //cout << "*****1" << endl;
            //cout << maxsize << endl;
            //cout << "*****1" << endl;
        }
        return;
    }

    if(flag == 0){
        //(A,Ca,B,Cb,0)
        int branchU = -1;
        int branchV = -1;


        vector<vertexDeg> degreeLIndex;
        vector<vertexDeg> vertexRDegIndex;
        bool ppf = true;
        while(ppf){
            auto it = newCa.begin();
            //for(auto &u:newCa){
            //cout<<"teszt1"<<endl;
            while(it!=newCa.end()){
                vector<int> *adju = &(bg.adjL[*it]);
                int degu = newB.size();
                //
                int size1 = adju->size();
                int size2 = newCb.size();
                if(size1<=size2){
                    for(auto &v: *adju){
                        if(newCb.find(v)!=newCb.end()){
                            degu++;
                        }
                    }
                }else{
                    for(auto &v:newCb){
                        if(bg.adjHL[*it].find(v)!=bg.adjHL[*it].end()){
                            degu++;
                        }
                    }
                }
                vertexDeg udegindex;
                udegindex.v = *it;
                udegindex.deg = degu;
                degreeLIndex.emplace_back(udegindex);
                //checkedU.emplace(u);
                bool f = true;
                bool f1 =true;
                if(degu<=maxsize){
                    //prunedU.emplace_back(u);
                    it = newCa.erase(it);
                    f=false;
                    f1=false;
                }

                //if((degu<(newB.size()+newCb.size()-2)) && degu>maxsize){
                //leftDegSat = false;
                //    branchU = u;
                //break;
                //}


                if(degu==newB.size()+newCb.size() && f1){
                    //optU.emplace(u);
                    newA.emplace(*it);
                    it = newCa.erase(it);
                    f=false;
                }

                //if(degu==newB.size()+newCb.size() && f1==false){
                //    cout<<"test potential problem";
                //}

                if(f){
                    it++;
                }

            }


            if(min((newA.size()+newCa.size()),(newB.size() + newCb.size()))<=maxsize){
                //cout<<"pruning 2"<<endl;
                return;
            }
            //after reduction, check if it can be terminated
            if(newCa.size() == 0 || newCb.size() == 0) {
                int currentsize = newA.size() + newCa.size();
                if (currentsize > newB.size() + newCb.size()) {
                    currentsize = newB.size() + newCb.size();
                }

                //cout<<end<<endl;
                if (currentsize > maxsize) {
                    maxsize = currentsize;
                    //cout << "*****2" << endl;
                    //cout << maxsize << endl;
                    //cout << "*****2" << endl;
                }
                return;
            }


            auto v = newCb.begin();
            bool f11 = true;

            //for (auto &v:newCb){
            while(v!=newCb.end()){
                vector<int> *adjv = &(bg.adjR[*v]);
                int degv = newA.size();

                int size1= adjv->size();
                int size2 = newCa.size();
                if(size1<= size2){
                    for(auto &u:*adjv){
                        if(newCa.find(u)!=newCa.end()){
                            degv++;
                        }
                    }
                }else{
                    for(auto &u:newCa){
                        if(bg.adjHR[*v].find(u)!=bg.adjHR[*v].end()){
                            degv++;
                        }
                    }

                }


                vertexDeg vdegind;
                vdegind.v = *v;
                vdegind.deg = degv;
                vertexRDegIndex.emplace_back(vdegind);
                bool f = true;
                bool f1 = true;

                if(degv<=maxsize){
                    //prunedV.emplace_back(v);
                    v=newCb.erase(v);
                    f = false;
                    f1= false;
                    ppf = true;
                    f11 = false;

                }
                //if((degv<(newA.size()+newCa.size()-2))&&degv>maxsize){
                //rightDegSat = false;
                //    branchV = v;
                // }
                if(degv==newA.size()+newCa.size()&& f1){
                    // optV.emplace(v);
                    newB.emplace(*v);
                    v=newCb.erase(v);
                    f = false;
                }

                //if(degv==newA.size()+newCa.size() && f1==false){
                //    cout<<"test potential problem 2";
                //}

                if(f){
                    v++;
                }

            }
            //if(ppf){

           // }
            if(f11==true){
                ppf = false;

            }else{
                //clean up degree index
                degreeLIndex.empty();
                vertexRDegIndex.empty();
                //cout<<"potential pruning 1"<<endl;
            }
            //cout<<"teszt4"<<endl;
            //for(auto &v:prunedV){
            //    newCb.erase(v);
            //}
            //for(auto &v:optV){
            //    newB.emplace(v);
            //    newCb.erase(v);
            // }
            //}


            if(min((newA.size()+newCa.size()),(newB.size() + newCb.size()))<=maxsize){
                //cout<<"pruning 3"<<endl;
                return;
            }


            if(newCa.size() == 0 || newCb.size() == 0) {
                int currentsize = newA.size() + newCa.size();
                if (currentsize > newB.size() + newCb.size()) {
                    currentsize = newB.size() + newCb.size();
                }

                //cout<<end<<endl;
                if (currentsize > maxsize) {
                    maxsize = currentsize;
                    //if(maxsize==10){
                    //    cout<<"for debug"<<endl;
                    // }
                    cout << "*****3" << endl;
                    cout << maxsize << endl;
                    cout << "*****3" << endl;
                }
                return;
            }

        }



        //compare degree with newA, newB, newCa, newCb again
        //vertexRDegIndex
        //degreeLIndex
        //branch at a vertex with currently minimum degree
        //computing degree in the subgraph cannot be pruned, which can also be used to compute the induced subgraph
        int minDegVec1;
        int minDegVec2 = -1;
        int minDegree = 100000;
        for(auto & udeg:degreeLIndex){
            if((udeg.deg<(newB.size()+newCb.size()-2))){
                //leftDegSat = false;
                branchU = udeg.v;
                if(udeg.deg<minDegree){
                    minDegree= udeg.deg;
                    minDegVec1 = branchU;
                }
                //break;
            }
        }



        for(auto & vcdeg:vertexRDegIndex){
            if((vcdeg.deg<(newA.size()+newCa.size()-2))){
                //rightDegSat = false;
                branchV = vcdeg.v;
                if(vcdeg.deg<minDegree){
                    minDegree = vcdeg.deg;
                    minDegVec2 = branchV;
                }
            }
        }



        if(branchV==-1&&branchU==-1){
            //p-time solvable

            unordered_map<int,unordered_set<int>> nadjL, nadjR, nadjLc,nadjRc;
            for(auto & u:newCa){
                unordered_set<int> adju;
                unordered_set<int> adjur;

                nadjL.emplace(u,adju);
                nadjLc.emplace(u,adjur);

            }

            for(auto &v:newCb){
                unordered_set<int> adjv;
                unordered_set<int> adjvr;

                nadjR.emplace(v,adjv);
                nadjRc.emplace(v,adjvr);

            }


            computeInG(nadjL,newCa,bg.adjHL, nadjR, newCb, bg.adjHR);
            computeInG(nadjR, newCb, bg.adjHR, nadjL,newCa,bg.adjHL);

            computeRInG(nadjLc, nadjL,newCb);
            computeRInG(nadjRc, nadjR,newCa);

            //cout<<"dense "<<endl;
            //printIG(nadjL,nadjR);
            //cout<<"sparse "<<endl;
            //printIG(nadjLc,nadjRc);



            //cout<<"dense "<<endl;
            //printIG(nadjR,nadjL);
            //cout<<"sparse "<<endl;
            //printIG(nadjRc,nadjLc);

            vector<pair<int,int>> paths;
            //cout<<"here! 2"<<endl;
            //void denseMBB::dynamicMBB(vector<pair<int,int>>& paths, unordered_set<int> & A,unordered_map<int, unordered_set<int>> &adjL, unordered_set<int> & B, unordered_map<int, unordered_set<int>> &adjR) {
            dynamicMBB(paths, newA,nadjLc,newB,nadjRc);




            return;
        }
        branchU = minDegVec1;
        branchV = minDegVec2;
        if(minDegVec2 == -1){
            //

            newA.emplace(branchU);
            newCa.erase(branchU);

            unordered_set<int> newnewCb;
            int sizeadj = bg.adjL[branchU].size();
            int sizeCb = newCb.size();
            if(sizeadj<=sizeCb){
                vector<int> *adju = &(bg.adjL[branchU]);
                for(auto & v:*adju){
                    if(newCb.find(v)!=newCb.end()){
                        newnewCb.emplace(v);
                    }
                }
            }else{
                for(auto & v:newCb){
                    if(bg.adjHL[branchU].find(v)!=bg.adjHL[branchU].end()){

                        newnewCb.emplace(v);
                    }
                }
            }


            recEnu(newB,newnewCb,newA,newCa,1);
            //recover pruned and opt
            newA.erase(branchU);

            recEnu(newA,newCa,newB,newCb,0);

            return;
        }else{
            newB.emplace(branchV);
            newCb.erase(branchV);

            unordered_set<int> newnewCa;

            int sizeadjv = bg.adjR[branchV].size();
            int sizeca = newCa.size();

            if(sizeadjv<=sizeca){

                for(auto &u:bg.adjR[branchV]){
                    if(newCa.find(u)!=newCa.end()){
                        newnewCa.emplace(u);
                    }
                }

            }else{
                for(auto &u:newCa){
                    if(bg.adjHR[branchV].find(u)!=bg.adjHR[branchV].end()){

                        newnewCa.emplace(u);
                    }
                }
            }

            recEnu(newA,newnewCa,newB,newCb,0);
            newB.erase(branchV);

            recEnu(newB,newCb,newA,newCa,1); //lead to 4
            return;
        }

    }else{
        //R
        //adjR
        //(A,Ca,B,Cb,1)
        // (B,Cb,A,Ca,1)
        //A is B in fact
        //Ca is Cb in fact


        //compute induced subgraph, check p-time solvable condition, and pruning
        //process R
        int branchV = -1;
        int branchU = -1;
        //vector<int> prunedV;
        //unordered_set<int>optV;
        vector<vertexDeg> cadegind;
        vector<vertexDeg> cbdegind;

        bool ppf = true;
        while(ppf){
            auto v = newCa.begin();
            //for(auto &v:newCa){
            //cout<<"teszt5"<<endl;



            while(v!=newCa.end()){
                vector<int> *adjv = &(bg.adjR[*v]);  // copy!! fixed
                int degv = newB.size();
                int size1 = adjv->size();
                int size2 = newCb.size();
                if(size1<=size2){
                    for(auto &u:*adjv){
                        if(newCb.find(u)!=newCb.end()){
                            degv++;
                        }
                    }

                }else{
                    for(auto &u:newCb){
                        if(bg.adjHR[*v].find(u)!=bg.adjHR[*v].end()){
                            degv++;
                        }
                    }
                }

                vertexDeg  vdeg;
                vdeg.v =*v ;
                vdeg.deg =  degv;
                cadegind.emplace_back(vdeg);

                bool f = true;
                bool f1 = true;
                if(degv<=maxsize){
                    //prunedV.emplace_back(*v);
                    v = newCa.erase(v);
                    f= false;
                    f1= false;
                }

                //if((degv<(newB.size()+newCb.size()-2))&&degv>maxsize){
                //rightDegSat = false;
                //    branchV = v;
                //}

                if(degv==newB.size()+newCb.size()&&f1){
                    //optV.emplace(v);
                    newA.emplace(*v);
                    v = newCa.erase(v);
                    f= false;
                }

                //if(degv==newB.size()+newCb.size() && (f1==false)){
                //    cout<<"potential problem 3"<<endl;
                //}

                if(f){
                    v++;
                }

            }
            //cout<<"teszt6"<<endl;
            //reduction for R

            //for(auto &v:prunedV){
            //   newCa.erase(v);
            //}

            //for(auto &v:optV){
            //    newA.emplace(v);
            //    newCa.erase(v);
            //}
            if(min((newA.size()+newCa.size()),(newB.size() + newCb.size()))<=maxsize){
                //cout<<"pruning 4"<<endl;
                return;
            }

            if(newCa.size() == 0 || newCb.size() == 0) {
                int currentsize = newA.size() + newCa.size();
                if (currentsize > newB.size() + newCb.size()) {
                    currentsize = newB.size() + newCb.size();
                }

                //cout<<end<<endl;
                if (currentsize > maxsize) {
                    maxsize = currentsize;
                    cout << "*****4" << endl;
                    cout << maxsize << endl;
                    cout << "*****4" << endl;
                }
                return;
            }

            //process L

            // vector<int> prunedU;
            // unordered_set<int>optU;

            auto u = newCb.begin();
            //cout<<"teszt7"<<endl;



            bool ppf2 = false;
            bool pp22 = false;
            while(u!=newCb.end()){
                //for(auto &u:newCb){
                vector<int> *adju = &(bg.adjL[*u]); //previous copy
                int degu = newA.size();
                int size1 = adju->size();
                int size2 = newCa.size();
                if(size1<=size2){
                    for(auto &v:*adju){
                        if(newCa.find(v)!=newCa.end()){
                            degu++;
                        }
                    }
                }else{
                    //cout<<", problem?"<<endl;
                    for(auto &v:newCa){
                        if(bg.adjHL[*u].find(v)!=bg.adjHL[*u].end()){
                            degu++;
                        }

                    }

                }

                vertexDeg udeg;
                udeg.v =*u;
                udeg.deg = degu;
                cbdegind.emplace_back(udeg);
                bool f = true;
                bool f1 = true;

                if(degu<=maxsize){
                    //prunedU.emplace_back(u);
                    u=newCb.erase(u);
                    f = false;
                    f1 = true;
                    ppf2 = true;
                    pp22 = true;
                }

                //if((degu<(newA.size()+newCa.size()-2)) && degu>maxsize){
                //    branchU = u;
                //}

                if(degu==newA.size()+newCa.size()&&f1){
                    //optU.emplace(*u);
                    newB.emplace(*u);
                    u=newCb.erase(u);
                    f = false;
                }
                //if(degu==newA.size()+newCa.size() && (f1==false)){
                //    cout<<"potential problem 4"<<endl;
                //}

                if(f){
                    u++;
                }

            }

            //if(ppf2){
            //    cout<<"potential pruning 2"<<endl;
            //}

            if(pp22 == false){
                ppf = false;
            }else{
                cadegind.empty();
                cbdegind.empty();

            }
        }


        //for(auto &u:prunedU){
        //    newCb.erase(u);
        //}

        //for(auto &u:optU){
        //    newB.emplace(u);
         //   newCb.erase(u);
       // }
        //cout<<"teszt8"<<endl;
        if(min((newA.size()+newCa.size()),(newB.size() + newCb.size()))<=maxsize){
            //cout<<"pruning 5"<<endl;
            return;
        }

        if(newCa.size() == 0 || newCb.size() == 0) {
            int currentsize = newA.size() + newCa.size();
            if (currentsize > newB.size() + newCb.size()) {
                currentsize = newB.size() + newCb.size();
            }

            //cout<<end<<endl;
            if (currentsize > maxsize) {
                maxsize = currentsize;
                cout << "*****5" << endl;
                cout << maxsize << endl;
                cout << "*****5" << endl;
            }
            return;
        }

        //cadegind
        int minDegVec1;
        int minDegVec2 = -1;
        int minDegree = 10000000;

        for(auto & udeg:cadegind){
            if((udeg.deg<(newB.size()+newCb.size()-2))){
                //rightDegSat = false;
                branchV = udeg.v;
                if(udeg.v<minDegree){
                    minDegree = udeg.v;
                    minDegVec1 = branchV;
                }
            }
        }

        //cbdegind


        for(auto & vdeg:cbdegind){
            if((vdeg.deg<(newA.size()+newCa.size()-2))){
                branchU = vdeg.v;
                if(vdeg.deg<minDegree){
                    minDegree = vdeg.deg;
                    minDegVec2 = branchU;
                }
            }

        }


        if(branchV==(-1)&&branchU==(-1)){
            //p-time solvable

            //unordered_map<int,unordered_set<int>> iadjL;
            //unordered_map<int,unordered_set<int>> iadjR;
            //cout<<"potential size: "<<min((newA.size()+newCa.size()),(newB.size() + newCb.size()))<<endl;
            //cout<<"p-time solvable!"<<endl;

            //newCa is in R in fact
            //newCb is in L in fact

            //
            unordered_map<int,unordered_set<int>> nadjL, nadjR, nadjLc, nadjRc;
            for(auto & v:newCa){
                unordered_set<int> adjv;
                unordered_set<int> adjvc;

                nadjL.emplace(v,adjv);
                nadjLc.emplace(v,adjvc);

            }

            for(auto &u:newCb){
                unordered_set<int> adju;
                unordered_set<int> adjuc;

                nadjR.emplace(u,adju);
                nadjRc.emplace(u,adjuc);
            }

            computeInG(nadjL,newCa,bg.adjHR, nadjR, newCb, bg.adjHL);
            computeInG(nadjR, newCb, bg.adjHL, nadjL,newCa,bg.adjHR);
            //computeInG(nadjL,newCa,bg.adjHL, nadjR, newCb, bg.adjHR);
            //computeInG(nadjR, newCb, bg.adjHR, nadjL,newCa,bg.adjHL);



            computeRInG(nadjLc, nadjL,newCb);
            computeRInG(nadjRc, nadjR,newCa);



            //cout<<"dense "<<endl;
           // printIG(nadjR,nadjL);
            //cout<<"sparse "<<endl;
            //printIG(nadjRc,nadjLc);

            vector<pair<int,int>> paths;
            cout<<"here! 2"<<endl;
            //void denseMBB::dynamicMBB(vector<pair<int,int>>& paths, unordered_set<int> & A,unordered_map<int, unordered_set<int>> &adjL, unordered_set<int> & B, unordered_map<int, unordered_set<int>> &adjR) {
            dynamicMBB(paths, newA,nadjLc,newB,nadjRc);
            //bfsLR(paths,nadjRc,nadjLc);


            return;
        }
        branchV = minDegVec1;
        branchU = minDegVec2;

        //reduction finish
        if(minDegVec2==-1){
            // (B,Cb,A,Ca,1)
            //(A,Ca,B,Cb,1)
            //branching at branchV
            newA.emplace(branchV);
            newCa.erase(branchV);

            unordered_set<int> newnewCb;

            int sizeadjv = bg.adjR[branchV].size();
            int sizeCa = newCa.size();

            if(sizeadjv<=sizeCa){
                for(auto &u:bg.adjR[branchV]){
                    if(newCa.find(u)!=newCa.end()){
                        newnewCb.emplace(u);
                    }
                }

            }else{
                for(auto &u:newCa){
                    if(bg.adjHR[branchV].find(u)!=bg.adjHR[branchV].end()){
                        newnewCb.emplace(u);
                    }
                }
            }
            //B in fact is A
            recEnu(newB,newnewCb,newA,newCa,0);
            newA.erase(branchV);


            recEnu(newA,newCa,newB,newCb,1);
            return;


        }else{
            // (B,Cb,A,Ca,1)
            //branch at branchU
            newB.emplace(branchU);
            newCb.erase(branchU);

            unordered_set<int> newnewCa;

            int sizeofu = bg.adjL[branchU].size();
            int sizeofCb = newCb.size();
            if(sizeofu<=sizeofCb){
                for(auto &v:bg.adjL[branchU]){
                    if(newCb.find(v) != newCb.end()){
                        newnewCa.emplace(v);
                    }
                }
            }else{
                for(auto &v:newCb){
                    if(bg.adjHL[branchU].find(v)!=bg.adjHL[branchU].end()){
                        newnewCa.emplace(v);
                    }
                }
            }

            recEnu(newA,newnewCa,newB,newCb,1); //lead to 4
            newB.erase(branchU);

            recEnu(newB,newCb,newA,newCa,0);
            return;
        }
    }



}

void denseMBB::indstart() {

    if(bg.adjL.size()<bg.adjR.size()){
        //sort vertex by degree
        //
        int maxdeg = 0;
        for(auto &ns:bg.adjL){
            if(ns.size()>maxdeg){
                maxdeg = ns.size();
            }
        }

        int *nd = new int [maxdeg+1]();

        for(auto &ns:bg.adjL){
            nd[ns.size()]=nd[ns.size()]+1;
        }

        int spos=0;
        for(int i =0; i<maxdeg+1; i++){
            int temp = nd[i];
            nd[i] = spos;
            spos = temp+spos;
        }

        int *sortedV = new int [bg.adjL.size()]();

        for(int i=0;i<bg.adjL.size();i++){
            sortedV[nd[bg.adjL[i].size()]]=i;
            nd[bg.adjL[i].size()]++;
        }
        delete [] nd;

       // for(int i=0;i<bg.adjL.size();i++){
        //    cout<<sortedV[i]<<"' degree is"<<bg.adjL[sortedV[i]].size()<<endl;
       // }


        delete [] sortedV;
    }else{

        int maxdeg = 0;
        for(auto &ns:bg.adjR){
            if(ns.size()>maxdeg){
                maxdeg = ns.size();
            }
        }

        int *nd = new int [maxdeg+1]();

        for(auto &ns:bg.adjR){
            nd[ns.size()]=nd[ns.size()]+1;
        }

        int spos=0;
        for(int i =0; i<maxdeg+1; i++){
            int temp = nd[i];
            nd[i] = spos;
            spos = temp+spos;
        }

        int *sortedV = new int [bg.adjR.size()]();

        for(int i=0;i<bg.adjR.size();i++){
            sortedV[nd[bg.adjR[i].size()]]=i;
            nd[bg.adjR[i].size()]++;
        }
        delete [] nd;

        //for(int i=0;i<bg.adjR.size();i++){
         //   cout<<sortedV[i]<<"' degree is"<<bg.adjR[sortedV[i]].size()<<endl;
        //}
        //we are dealing with B
        int count = 0;
        for(int i=bg.adjR.size()-1; i>=0;i--){
        //for(int i=0;i<bg.adjR.size();i++){
            count++;
            int branchingV = sortedV[i];
            cout<<"current processing v is "<< branchingV<<". It is "<<count<<"th vertices out of "<<bg.adjR.size()<<endl;
            unordered_set<int> B;
            unordered_set<int> Cb;

            unordered_set<int> A;
            unordered_set<int> Ca;



            //R
            B.emplace(branchingV);

            for(int j=i-1;j>=0;j--){
            //for(int j=i+1;j<bg.adjR.size();j++){
                Cb.emplace(sortedV[j]);
            }

            for(auto& u: bg.adjR[branchingV]){
                Ca.emplace(u);
            }

            unordered_set<int> newCb;

            for(auto &u:Ca){
                for(auto &v:bg.adjL[u]){
                    if(v!=branchingV){
                        newCb.emplace(v);
                    }

                }
            }

            cout<<"newCb size "<<newCb.size() <<endl;
            cout<<"cb size "<<Cb.size() <<endl;
            //cout<<"alg starts degree"<<endl;
            //if(Ca.size()>maxsize){
            recEnu(A,Ca,B,Cb,0);
            //}

        }
        delete [] sortedV;
    }

}


//all use relative logic
void denseMBB::computeInG(unordered_map<int, unordered_set<int>>  &nadjL, unordered_set<int>& L, vector<unordered_set<int>> &adjL,
                          unordered_map<int, unordered_set<int>>  &nadjR, unordered_set<int>& R, vector<unordered_set<int>> &adjR) {

    //nadjL must be a subset of R
    for(auto u:L){
        unordered_set<int> &adju = adjL[u];
        if(adju.size()<=R.size()){
            for(auto v:adju){
                if(R.find(v)!=R.end()){
                    nadjL[u].emplace(v);
                }
            }
        }else{
            for(auto v:R){
                if(adju.find(v)!=adju.end()){
                    nadjL[u].emplace(v);
                }
            }

        }
    }

}



void denseMBB::printIG(unordered_map<int, unordered_set<int>> & L, unordered_map<int, unordered_set<int>> & R) {
    cout<<"adj for L is: "<<endl;
    for(auto v:L){
        cout<<v.first<<":";
        for(auto u:v.second){
            cout<<u<<",";
        }
        cout<<endl;
    }


    cout<<"adj for R is: "<<endl;
    for(auto v:R){
        cout<<v.first<<":";
        for(auto u:v.second){
            cout<<u<<",";
        }
        cout<<endl;
    }

}

void denseMBB::computeRInG(unordered_map<int,unordered_set<int>> & nadjL, unordered_map<int, unordered_set<int>> & adjL, unordered_set<int> & R) {

    for(auto &u: adjL){
        for(auto &v:R){
            if(u.second.find(v)==u.second.end()){
                nadjL[u.first].emplace(v);
            }
        }
    }


}

void denseMBB::dynamicMBB(vector<pair<int,int>>& paths, unordered_set<int> & A,unordered_map<int, unordered_set<int>> &adjL, unordered_set<int> & B, unordered_map<int, unordered_set<int>> &adjR) {
    //collect paths or cycle
    // adjL and adjR represent the complementary subgraph
    bfsLR(paths,adjL,adjR); //get the paths

    //now paths are here
    //a pair in the paths

    //generate sizes
    vector<vector<pair<int,int>>> sizesForP;
    //
    vector<pair<int,int>> sizes;
    pair<int,int> size0(A.size(),B.size());
    sizes.push_back(size0);
    sizesForP.emplace_back(sizes);
    for(auto &p:paths){
        vector<pair<int,int>> sizes;
        calcuateSizes(p,sizes); //p.first is the type, p.second is the length
        sizesForP.emplace_back(sizes);
    }

    int sizeL = A.size()+adjL.size();
    int sizeR = B.size()+adjR.size();

    int** t = new int*[sizeL+1];
    vector<vector<pair<int,int>>> tprime;
    for (int i = 0; i <=sizeL; i++){
        t[i] = new int[sizeR+1]{0};
        vector<pair<int,int>> a(sizeR+1);
        tprime.emplace_back(a);
    }


    t[A.size()][B.size()] = 1;
    tprime[A.size()][B.size()].first = A.size();
    tprime[A.size()][B.size()].second = B.size();
    vector<pair<int,int> *> intermediateSizes;
    intermediateSizes.emplace_back(&(tprime[A.size()][B.size()]));
    vector<pair<int,int> *> *tt = &intermediateSizes;

    //argmax
    for(int i=1;i<sizesForP.size();i++){
        vector<pair<int,int> *> newintermediateSizes;
        for(auto &size:*tt){
            for(auto &sizei:sizesForP[i]){
                if(t[size->first+sizei.first][size->second+sizei.second]!=i+1){
                    t[size->first+sizei.first][size->second+sizei.second] = i+1;
                    int currentSize = min(size->first+sizei.first,size->second+sizei.second);
                    if(currentSize>maxsize){
                            maxsize = currentSize;
                    }
                    pair<int,int>ab(sizei.first,sizei.second);
                    tprime[size->first+sizei.first][size->second+sizei.second] = ab;
                    newintermediateSizes.emplace_back(&ab);
                }
            }
        }
        *tt = newintermediateSizes; //now it is copy, shall be replaced by a pointer. However, it would not affect the time complexity
    }

    //print t

    //cout<<"max size is "<<maxsize<<"************************"<<endl;
    //design your own biclique collection me

    delete [] t;
}

void denseMBB::bfsLR(vector<pair<int,int>>&paths,unordered_map<int, unordered_set<int>> &adjL, unordered_map<int, unordered_set<int>> &adjR) {
    unordered_set<int> visitedL, visitedR;
    int flag = 0; // 0->left, 1->right
    stack<pair<int,sInf>> stack;

    for(auto v:adjL){ //start from L
        if(visitedL.find(v.first)==visitedL.end()){
            int minDegree = adjL[v.first].size(); //the minimum degree of the current path
            //put it into stack
            pair<int,sInf> s;
            sInf si;
            si.flag = 0;
            si.vistis = 0;
            s = make_pair(v.first,si); // 0 means v is in the L
            stack.push(s); // s is copied, so s can be reused later
            int length=0;
            while(!stack.empty()){
                s = stack.top();
                stack.pop();
                if(s.second.flag==0){
                    // s.first is a vertex in L
                    if(visitedL.find(s.first)==visitedL.end()){
                    //if( s.second.vistis==0){

                        //length++;
                        visitedL.emplace(s.first);
                        for(auto u: adjL[s.first]){ // u is in R

                           if(visitedR.find(u)==visitedR.end()){
                               if(adjR[u].size()<minDegree){
                                   minDegree = adjR[u].size();
                               }
                               length++;
                           }
                        }
                    }

                }

                if(s.second.flag==1){
                    // s.first is a vertex in R
                    if(visitedR.find(s.first)==visitedR.end()){
                        visitedR.emplace(s.first);

                        for(auto v: adjR[s.first]){ // v is in L
                            if(visitedL.find(v)==visitedL.end()){
                                if(adjL[v].size()<minDegree){
                                    minDegree = adjL[v].size();
                                }
                                length++;
                            }

                        }

                    }

                }

            }
            if(minDegree==1){
                //it is a path
                cout<<"a path with length of "<<length<<endl;
                pair<int,int> s;
                if(length%2==0){//even path
                    // two cases for an even path
                    if(adjL.size()>adjR.size()){
                        s.first=2;
                        s.second=length;
                    }else{
                        s.first=3;
                        s.second=length;
                    }

                }else{//odd path
                    s.first=1;
                    s.second=length;
                }

                paths.push_back(s);
            }

            if(minDegree ==2){
                //it is a cycle
                cout<<"a cycle with length of "<<length<<endl;
                pair<int,int> s(0,length);
                paths.push_back(s);
            }

        }


    }

}

void denseMBB::calcuateSizes(pair<int, int> & p, vector<pair<int,int>> & sizes) {

    if(p.first ==0){  //cycle
        div_t result;
        pair<int,int> size1, size2;
        result = div(p.second,2);
        size1.first = 0;
        size1.second = result.quot;

        size2.first = result.quot;
        size2.second = 0;

        if(p.second==4){
            sizes.emplace_back(size1);
            sizes.emplace_back(size2);

        }else{
            result = div(p.second,2);
            int divR= result.quot;
            int loopEnd = divR - 2;

            //loop for the invariance
            sizes.emplace_back(size1);
            for(int i=1; i<=loopEnd; i++){
                pair<int,int> size;
                size.first = i;
                size.second = divR-i-1;
                sizes.emplace_back(size);
            }
            sizes.emplace_back(size2);
        }

    }

    if(p.first==1){
        //odd path
        div_t result;
        result = div(p.second+1,2);
        int divR = result.quot;
        for(int i=0;i<=divR;i++){
            pair<int,int> size;
            size.first = i;
            size.second = divR - i;
            sizes.emplace_back(size);
            //cout<<"size is: "<< size.first<<", "<<size.second<<endl;
        }
        //cout<<"done!"<<endl;
    }

    if(p.first==2){
        //even path L>R
        div_t result;
        result = div(p.second,2);
        int divR= result.quot;
        for(int i=0; i<=divR-1; i++){
            pair<int,int> size;
            size.first = i;
            size.second = divR-i;
            sizes.emplace_back(size);
        }
        pair<int,int> sizeLast;
        sizeLast.first = divR+1;
        sizeLast.second = 0;
        sizes.emplace_back(sizeLast);
    }

    if(p.first==3){
        //even path L<R
        div_t result;
        result = div(p.second,2);
        int divR= result.quot;
        for(int i=0; i<=divR-1; i++){
            pair<int,int> size;
            size.first = divR-i;
            size.second =  i;
            sizes.emplace_back(size);
        }
        pair<int,int> sizeLast;
        sizeLast.first = 0 ;
        sizeLast.second = divR+1;
        sizes.emplace_back(sizeLast);
    }


}





