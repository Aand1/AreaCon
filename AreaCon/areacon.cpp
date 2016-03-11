/********************************************/
/**
* @mainpage AreaCon: A C++ Library for Area-Constrained Partitioning
* @details AreaCon is a light-weight, C++ library for carrying out area-constrained partitioning operations in floating-point precision. The library is primarily a numerical implementation of the area-constrained partitioning algorithm by Patel, Frasca, and Bullo, 2014 (see http://www.areacon.org/References ).

Low-level source code documentation is provided here; however, additional discussion about how to use the code along with numerical examples can be found on the project's Wiki ( https://github.com/jrpeters/AreaCon/wiki ).

Please cite AreaCon whenever possible. A bibtex citation is provided below:
@verbatim
 @Misc{AreaCon:16,
 author =        {J.R. Peters and Contributors},
 title =         {The {AreaCon} Library},
 howpublished =  {\texttt{http://www.areacon.org}},
 year =          {2016},
 note =          {v. 1},
 abstract =      {A C++ library for area-constrained partitioning.},
 }
@endverbatim

For more information on the code and information on how to contact me, please see http://www.areacon.org
* @author Jeffrey R. Peters ( http://www.jeffreyrpeters.com , https://github.com/jrpeters )
* @version 1.0
* @date Feb. 2016
* @copyright Copyright &copy; 2016. The Regents of the University of California. All rights reserved. Licensed pursuant to the terms and conditions available for viewing at: http://opensource.org/licenses/BSD-3-Clause .

* @copyright The Clipper library ( www.angusj.com/delphi/clipper.php ), as used by AreaCon, is also subject to the following:

* @copyright Copyright &copy; 2016. The Regents of the University of California. Distributed under the Boost Software License, Version 1.0. (See http://www.boost.org/LICENSE_1_0.txt ).
***********************************************/
#include "areacon.h"

namespace AreaCon {
    // Point Class-------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    double Point::Robustness_Constant;
    Point::Point(const double x, const double y):x(x), y(y){}
    double Point::Norm() const {return Norm(*this);}
    Point Point::AddPoint(const Point Test) const {return AddPoints(*this, Test);}
    Point Point::FindPerpDirection(const Point Test, const double Norm) const {return FindPerpDirection(*this, Test, Norm);}
    void Point::FlipDirection(void){x = -x; y = -y;}
    double Point::PerpDistanceToLine(const Point Test1, const Point Test2) const{return PerpDistanceToLine(Test1, Test2, *this);}
    bool Point::AreCollinear(const Point Test1, const Point Test2)const{ return AreCollinear(Test1,Test2, *this);}
    bool Point::AreBetween(const Point Test1, const Point Test2) const{return AreBetween(Test1,Test2,*this);}
    void Point::Mult(const double factor){x = factor*x, y =  factor*y;}
    bool Point::IsEqual(const Point Test1, const Point Test2){if(Test1.x == Test2.x && Test1.y == Test2.y){return true;}else{return false;}}
    double Point::Distance(const Point Test1,const Point Test2) {return sqrt((Test1.x-Test2.x)*(Test1.x-Test2.x)+(Test1.y-Test2.y)*(Test1.y-Test2.y));}
    Point Point::FindPointAlongLine(const Point Test1, const Point Test2, const double distance){return Point(Test1.x+(Test2.x-Test1.x)*distance, Test1.y+(Test2.y-Test1.y)*distance);}
    double Point::Norm(const Point Test){return sqrt(Test.x*Test.x+Test.y*Test.y);}
    Point Point::AddPoints(const Point Test1, const Point Test2){
        return Point(Test1.x+Test2.x, Test1.y+Test2.y);
    }
    Point Point::FindPerpDirection(const Point Test1, const Point Test2, const double Norm){
        double distance = Distance(Test1, Test2);
        if (distance == 0){
            return Point(0,0);
        }else{
            return Point((Test2.y-Test1.y)/distance*Norm, (Test1.x-Test2.x)/distance*Norm);
        }
    }
    double Point::PerpDistanceToLine(const Point Test1, const Point Test2, const Point Test3){
        if (fabs(Test2.y-Test1.y)<Robustness_Constant){
            return fabs(Test3.y-Test2.y);
        }else if(fabs(Test2.x-Test1.x)<Robustness_Constant){
            return fabs(Test3.x-Test2.x);
        }else{
            return fabs((Test2.y-Test1.y)*Test3.x-(Test2.x-Test1.x)*Test3.y+Test2.x*Test1.y-Test2.y*Test1.x)/Distance(Test1, Test2);
        }
        
    }
    bool Point::AreCollinear(const Point Test1, const Point Test2, const Point Test3){
        double max_dist = Distance(Test1, Test2), temp;
        const double tolerance = Robustness_Constant;
        temp =Distance(Test2, Test3);
        if (temp>max_dist){
            max_dist = temp;
        }
        temp =Distance(Test1, Test3);
        if (temp>max_dist){
            max_dist = temp;
        }
        if (Test3.PerpDistanceToLine(Test1, Test2)/max_dist<tolerance){
            return true;
        }else{
            return false;
        }
    }
    bool Point::AreBetween(const Point Test1, const Point Test2, const Point Test3){
        if (!AreCollinear(Test1, Test2, Test3)){
            return false;
        }else{
            double distance = Distance(Test3, Test1)/Distance(Test2,Test1);
            if (distance>1 || distance <0){
                return false;
            }else if ( (Test3.x-Test1.x<=0 && Test2.x-Test1.x<=0)|| (Test3.x-Test1.x>=0 && Test2.x-Test1.x>=0)){
                return true;
                
            }else{
                return false;
            }
        }
    }
    std::vector<Point> Point::FindCollinearIntersection(const Point p1, const Point p2, const Point p3, const Point p4){
        std::vector<Point>  result;
        double tolerance = Robustness_Constant;
        if (AreBetween(p1,p2,p3)){
            result.push_back(p3);
            if (AreBetween(p1,p2,p4)){
                result.push_back(p4);
            }else if (AreBetween(p3, p4, p1) && Distance(p3, p1)>tolerance){
                result.push_back(p1);
            }else if (AreBetween(p3,p4,p2)&&Distance(p3,p2)>tolerance){
                result.push_back(p2);
            }
        }else if (AreBetween(p1,p2,p4)){
            result.push_back(p4);
            if (AreBetween(p3,p4,p1)&& Distance(p4,p1) >tolerance){
                result.push_back(p1);
            }else if (AreBetween(p3,p4,p2)&& Distance(p2, p4)>tolerance){
                result.push_back(p2);
            }
            
        }else if (AreBetween(p3,p4,p1) && AreBetween(p3,p4,p2)){
            result.push_back(p1);
            result.push_back(p2);
        }
        return result;
    }
    
// Poly Class-------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------
    
    Poly::Poly(std::vector<Point> Vertices):Vertices(Vertices), NPoly((int) Vertices.size()){InitializePoly();}
    std::vector<Point> Poly::GetVertices(void) const{return Vertices;}
    int Poly::GetNVertices(void) const{return NPoly;}
    void Poly::GetExtrema(double &minx, double &miny, double &maxx, double &maxy) const{minx = this->minx;maxx = this->maxx;miny = this->miny;maxy = this->maxy;}
    void Poly::SetVertices(const std::vector<Point> Vertices, const bool GetExtrema){
        this->Vertices = Vertices;
        this->NPoly = (int) Vertices.size();
        if (GetExtrema){
            InitializePoly();
        }else{
            bool flag = Vertices.empty();
            if (!flag && NPoly<3){
                throw std::runtime_error("List of vertices must contain at least 3 points");
            }
        }
    }
    bool Poly::pnpoly(const Point Test) const{
        bool value = false;
        if (Vertices.empty()){
            throw std::runtime_error("Polygon vertices have not been initialized");
        }else{
            Point Test_1 = Vertices.back(), Test_2;
            for (int ii = 0; ii<NPoly; ii++){
                Test_2 = Vertices[ii];
                if (Point::AreBetween(Test_1, Test_2, Test)){
                    value = true;
                    break;
                }
                if ((( Test_1.y <  Test.y && Test.y <= Test_2.y  ) ||( Test.y <= Test_1.y   && Test_2.y   < Test.y )) && ( Test_1.x   <=  Test.x || Test_2.x <= Test.x  )){
                    if ( Test_1.x  + ( Test.y - Test_1.y ) * ( Test_2.x - Test_1.x ) / ( Test_2.y - Test_1.y ) < Test.x )
                    {
                        value = !value;
                    }
                }
                Test_1 = Test_2;
                
            }
        }
        return value;
    }
    void Poly::InitializePoly(void){
        minx = INFINITY;
        maxx = -INFINITY;
        miny = minx;
        maxy = maxx;
        bool flag = Vertices.empty();
        if (!flag && NPoly<3){
            throw std::runtime_error("List of vertices must contain at least 3 points");
        }else if(!flag){
            for (int ii = 0; ii<NPoly; ii++){
                if (isinf(Vertices[ii].x) || isinf(Vertices[ii].y)){
                    throw std::runtime_error("Polygon vertices cannot be infinite");
                }else{
                    if (Vertices[ii].x<minx){
                        minx = Vertices[ii].x;
                    }
                    if (Vertices[ii].x>maxx){
                        maxx = Vertices[ii].x;
                    }
                    if (Vertices[ii].y<miny){
                        miny = Vertices[ii].y;
                    }
                    if (Vertices[ii].y>maxy){
                        maxy = Vertices[ii].y;
                    }
                    
                    for (int jj = ii+1; jj<NPoly; jj++){
                        if (Point::IsEqual(Vertices[ii],Vertices[jj])){
                            throw std::runtime_error("Polygon Vertices must all be distinct");
                        }
                    }
                }
            }
            if (minx == maxx ||miny == maxy){
                throw std::runtime_error("Polygon must have non-zero nominal area");
            }
        }
    }
    
    // Mult_Array Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------

    Mult_Array::Mult_Array(const int N):N(N){
        Array = new double*[N];
        for (int ii = 0; ii<N; ii++){
            Array[ii] = new double[N];
        }
    }
    Mult_Array::Mult_Array(const Mult_Array &obj):N(obj.N){
        Array = new double* [N];
        for (int ii = 0; ii<N; ii++){
            Array[ii] = new double [N];
            for (int jj = 0; jj<N; jj++){
                Array[ii][jj] = obj.Array[ii][jj];
            }
        }
    }
    Mult_Array& Mult_Array::operator=(const Mult_Array &obj){
        if (this != &obj){
            if (obj.N !=this->N){
                throw std::runtime_error("Incompatible sizes");
            }
            for (int ii = 0; ii<this->N; ii++){
                for (int jj = 0; jj<this->N;jj++){
                    this->Array[ii][jj] = obj.Array[ii][jj];
                }
            }
        }
        return *this;
    }
    Mult_Array::~Mult_Array(){
        for (int ii = 0; ii < N; ++ii){
            delete [] Array[ii];
        }
        delete [] Array;
    }
    // Delaunay Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    
    DelaunayGraph::DelaunayGraph(int NRegions):NRegions(NRegions){
        Graph = new Point** [NRegions];
        for (int ii = 0; ii<NRegions; ii++){
            Graph[ii] = new Point* [NRegions];
            for (int jj = 0;jj<NRegions;jj++){
                Graph[ii][jj] = new Point[2];
            }
        }
    }
    DelaunayGraph::~DelaunayGraph(){
        int N = NRegions;
        for (int ii = 0; ii < N; ++ii){
            for (int jj = 0; jj<N;++jj){
                delete [] Graph[ii][jj];
            }
            delete [] Graph[ii];
        }
        delete [] Graph;
    }
    DelaunayGraph::DelaunayGraph(const DelaunayGraph &obj):NRegions(obj.NRegions){
        Graph = new Point** [NRegions];
        for (int ii = 0; ii<NRegions; ii++){
            Graph[ii] = new Point* [NRegions];
            for (int jj = 0;jj<NRegions;jj++){
                Graph[ii][jj] = new Point[2];
                Graph[ii][jj][0] = obj.Graph[ii][jj][0];
                Graph[ii][jj][1] = obj.Graph[ii][jj][1];
            }
        }
    }
    DelaunayGraph& DelaunayGraph::operator=(const DelaunayGraph &obj){
        if (this != &obj) // protect against invalid self-assignment
        {
            if (this->NRegions != obj.NRegions){
                throw std::runtime_error("Incompatible Dimensions");
            }
            for (int ii = 0; ii<NRegions; ii++){
                for (int jj = 0;jj<NRegions;jj++){
                    this->Graph[ii][jj][0] = obj.Graph[ii][jj][0];
                    this->Graph[ii][jj][1] = obj.Graph[ii][jj][1];
                }
            }
        }
        
        return *this;
    }
    // Int_Params Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    
    Int_Params::Int_Params(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, std::vector<double> Int, std::vector<double> Intx, std::vector<double> Inty, double Unweighted_Area ):Coefficient_a(a), Coefficient_b(b), Coefficient_c(c), Coefficient_d(d), Int(Int), Intx(Intx), Inty(Inty), Unweighted_Area(Unweighted_Area){
        CheckParameters();
    }
    void Int_Params::CheckParameters(void) const{
        int size = (int) Coefficient_a.size();
        int size2 = (int) Int.size();
        
        if (size != Coefficient_b.size()){
            throw std::runtime_error("All integral coefficient vectors must be the same size");
        }else if (size != Coefficient_c.size()){
            throw std::runtime_error("All integral coefficient vectors must be the same size");
        }else if (size != Coefficient_d.size()){
            throw std::runtime_error("All integral coefficient vectors must be the same size");
        }else if (size2 !=0 && size2 !=size){
            throw std::runtime_error("If itegral vector is specified, it must be the same size as the integral coefficient vectors");
        }else if (size2 != Intx.size()){
            throw std::runtime_error("All integral vectors must be the same size");
        }else if (size2 != Inty.size()){
            throw std::runtime_error("All integral vectors must be the same size");
        }
    }
    // Parameters Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    Parameters::Parameters(const double line_int_step,const double weights_step, const double centers_step,const double volume_tolerance, const double convergence_criterion, const int max_iterations_volume, const int max_iterations_centers, const double Volume_Lower_Bound, const double Robustness_Constant):line_int_step(line_int_step), weights_step(weights_step), centers_step(centers_step), volume_tolerance(volume_tolerance), convergence_criterion(convergence_criterion), max_iterations_volume(max_iterations_volume), max_iterations_centers(max_iterations_centers), Volume_Lower_Bound(Volume_Lower_Bound), Robustness_Constant(Robustness_Constant){CheckParameters();};
    void Parameters::CheckParameters(void){
        if (line_int_step<=0){
            throw std::runtime_error("line_int_step must be greater than 0");
        }else if (weights_step<=0){
            throw std::runtime_error("weights_step must be greater than 0");
        }else if (centers_step<=0 || centers_step>1){
            throw std::runtime_error("centers_step must be greater than 0 and less than or equal to 1");
        }else if (volume_tolerance<=0){
            throw std::runtime_error("volume_tolerance must be greater than 0");
        }else if (convergence_criterion<=0){
            throw std::runtime_error("convergence_criterion must be greater than 0");
        }else if (max_iterations_volume<=0){
            throw std::runtime_error("max_iterations_volume must be greater than 0");
        }else if (max_iterations_centers<=0){
            throw std::runtime_error("max_iterations_centers must be greater than 0");
        }else if(Volume_Lower_Bound <= 0 || Volume_Lower_Bound>=1){
            throw std::runtime_error("Volume_Lower_Bound must be between 0 and 1");
        }else if (Robustness_Constant<=0){
            throw std::runtime_error("Robustness_Constant must be greater than 0");
        }
    }
    // Density Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    Density::Density():Volume_Lower_Bound(0){SetNewRegion(Region);}
    Density::Density(const Poly Region, const int Nx, const int Ny, const std::vector<double> Values):Volume_Lower_Bound(0){SetNewRegion(Region,Nx,Ny,Values);}
    void Density::SetExtrema(){Region.GetExtrema(minx, miny, maxx, maxy);}
    void Density::SetNewRegion(const Poly Region, const int Nx,const int Ny, const std::vector<double> Values){this->Region = Region;SetExtrema();SetParameters(Nx,Ny,Values);}
    Point Density::ConvertIndextoWorld(const int ii) const{int Index_x = ii/Ny, Index_y=ii%Ny;return Point(minx+Index_x*dx, miny+Index_y*dy);}
    void Density::CheckParameterSizes(void){
        int npoly = Region.GetNVertices();
        if (npoly == 0){
            Nx = 0;
            Ny = 0;
            Values = {};
        }else if (Nx == 0 || Ny == 0){
            Nx = 0;
            Ny = 0;
            Values = {};
        }else if (Nx*Ny != Values.size()){
            throw std::runtime_error("The size of Values must be equal to Nx*Ny");
        }
    }
    
    void Density::SetVolumeLowerBound(const double VolumeLowerBound){
        this->Volume_Lower_Bound = VolumeLowerBound;
    }
    double Density::GetVolumeLowerBound(void){
        return Volume_Lower_Bound;
    }
    
    
    void Density::Setdxy(void){
        double minx,miny,maxx,maxy;
        Region.GetExtrema(minx, miny, maxx, maxy);
        if (Nx!=0 && Ny!=0){
            dx = (maxx-minx)/(Nx-1);
            dy = (maxy-miny)/(Ny-1);
        }else{
            dx = 0;
            dy = 0;
        }
    }
    void Density::SetParameters(const int Nx, const int Ny, const std::vector<double> Values){
        this -> Nx = Nx;
        this -> Ny = Ny;
        this -> Values = Values;
        CheckParameterSizes();
        Setdxy();
        if (!Values.empty()){
            PreprocessIntegral();
        }
    }
    double Density::InterpolateValue(const Point &Test) const{
        int i = (int) ((Test.x-minx)/dx), j = (int) ((Test.y-miny)/dy);
        int i1,j1;
        Point p = ConvertIndextoWorld(i*Ny+j);
        double val00, val01, val10, val11;
        double val0s,valr0,valr1,val1s;
        double ys,xr;
        
        if (i == Nx-1){
            xr = 0;
            i1 = i;
        }else{
            xr = (Test.x-p.x)/dx;
            i1 = i+1;
        }
        if (j == Ny-1){
            ys = 0;
            j1 = j;
        }else{
            ys = (Test.y-p.y)/dy;
            j1 = j+1;
        }
        
        val00 = Values[Ny*i+j];
        val10 = Values[Ny*i1+j];
        val01 = Values[Ny*i+j1];
        val11 = Values[Ny*i1+j1];
        
        val0s = val00+(val01-val00)*ys;
        valr1 = val01+(val11-val01)*xr;
        val1s = val10+(val11-val10)*ys;
        valr0 = val00+(val10-val00)*xr;
        
        return (val0s+valr0-val00)+xr*(val00+val1s-val0s-val10)+ys*(val00+valr0-val01-valr1)+xr*ys*(val01+val10-val00-val11);
    }
    void Density::CreateIntegralCoefficients(void){
        Int_Params New;
        Integral = New;
        GridInRegion.clear();
        
        double a = 0, b = 0, c = 0, d = 0, gamma = 0, eta = 0, xi = 0,yval = miny, xval = minx;
        for (int ii = 0; ii< Nx-1; ii++){
            yval = miny;
            for (int jj = 0; jj<Ny-1;jj++){
                GridInRegion.push_back(Region.pnpoly(ConvertIndextoWorld(ii*Ny+jj)));
                gamma = -(1/(dx*dy))*(Values[ii*Ny+jj]+Values[(ii+1)*Ny+jj]-Values[ii*Ny+jj]-Values[(ii+1)*Ny+jj+1]);
                eta = (Values[(ii+1)*Ny+jj]-Values[ii*Ny+jj])/dx;
                xi = -(Values[ii*Ny+jj]-Values[ii*Ny+jj+1])/dy;
                a = -gamma*yval+eta;
                b = -gamma*xval+xi;
                c = gamma;
                d = xval*yval*gamma-yval*xi-xval*eta+Values[ii*Ny+jj];
                Integral.Coefficient_a.push_back(a);
                Integral.Coefficient_b.push_back(b);
                Integral.Coefficient_c.push_back(c);
                Integral.Coefficient_d.push_back(d);
                
                if (jj == Ny-2){
                    GridInRegion.push_back(Region.pnpoly(ConvertIndextoWorld(ii*Ny+jj+1)));
                }else{
                    yval += dy;
                }
            }
            if (ii == Nx-2){
                yval = miny;
                for (int jj = 0;jj<Ny;jj++){
                    GridInRegion.push_back(Region.pnpoly(ConvertIndextoWorld((ii+1)*Ny+jj)));
                    yval += dy;
                }
            }else{
                xval += dx;
            }
            
            
        }
    }
    double Density::CreateIntegralVector(void){
        Integral.Unweighted_Area = 0;
        double result = 0, resultx = 0,resulty = 0, xval,xval1,yval,yval1, total=0;
        int index = 0;
        xval1 = minx;
        for (int ii = 0;ii<Nx-1;ii++){
            xval = xval1;
            xval1 = xval+dx;
            yval1 = miny;
            for (int jj = 0;jj<Ny-1;jj++){
                yval = yval1;
                yval1 = yval+dy;
                result = dy*dx*Integral.Coefficient_d[index]+dy*(xval1*xval1-xval*xval)/2*Integral.Coefficient_a[index]+dx*(yval1*yval1-yval*yval)/2*Integral.Coefficient_b[index]+(yval1*yval1-yval*yval)*(xval1*xval1-xval*xval)/4*Integral.Coefficient_c[index];
                resultx = dy*(xval1*xval1-xval*xval)*Integral.Coefficient_d[index]/2+dy*(xval1*xval1*xval1-xval*xval*xval)*Integral.Coefficient_a[index]/3+(xval1*xval1-xval*xval)*(yval1*yval1-yval*yval)/4*Integral.Coefficient_b[index]+(yval1*yval1-yval*yval)*(xval1*xval1*xval1-xval*xval*xval)/6*Integral.Coefficient_c[index];
                resulty = (yval1*yval1-yval*yval)*dx*Integral.Coefficient_d[index]/2+(yval1*yval1-yval*yval)*(xval1*xval1-xval*xval)/4*Integral.Coefficient_a[index]+dx*(yval1*yval1*yval1-yval*yval*yval)/3*Integral.Coefficient_b[index]+(yval1*yval1*yval1-yval*yval*yval)*(xval1*xval1-xval*xval)/6*Integral.Coefficient_c[index];
                Integral.Int.push_back(result);
                Integral.Intx.push_back(resultx);
                Integral.Inty.push_back(resulty);
                if (Region.pnpoly(ConvertIndextoWorld(ii*Ny+jj))){
                    if (Region.pnpoly(ConvertIndextoWorld((ii+1)*Ny+jj))){
                        if (Region.pnpoly(ConvertIndextoWorld(ii*Ny+jj+1))){
                            if (Region.pnpoly(ConvertIndextoWorld((ii+1)*Ny+jj+1))){
                                total+= result;
                                Integral.Unweighted_Area +=dx*dy;
                            }
                        }
                    }
                }
                index++;
            }
        }
        return total;
    }
    void Density::NormalizeIntegralVector(const double &Total){
        if (Total == 0){
            std::cout<<"Warning: Density values do not have sufficient support. Treated as Uniform"<<std::endl;
            for (int ii = 0; ii<(Nx)*(Ny); ii++){
                Values[ii] = 1/Integral.Unweighted_Area;
            }
            SetParameters(Nx, Ny, Values);
        }else{
            for (int ii = 0; ii<(Nx-1)*(Ny-1); ii++){
                Integral.Int[ii] /= Total;
                Integral.Intx[ii] /= Total;
                Integral.Inty[ii] /= Total;
            }
        }
        
        
        
    }
    void Density::PreprocessIntegral(void){
        double Total;
        CreateIntegralCoefficients();
        Total = CreateIntegralVector();
        NormalizeIntegralVector(Total);
    }
    double Density::LineIntegral(double spacing, const Point &p1, const Point &p2) const{
        if (Values.empty()){
            throw std::runtime_error("Values have not been set!");
        }else if(spacing <= 0 || spacing >1){
            throw std::runtime_error("Spacing cannot be less than or equal to 0 or greater than 1");
        }
        double sum = 0, previous = 0, previous2 = 0;
        Point Test;
        previous = InterpolateValue(p1);
        for (double ii = 0;ii<1-spacing;ii+=spacing){
            Test.x = p1.x+(p2.x-p1.x)*(ii+spacing);
            Test.y = p1.y + (p2.y-p1.y)*(ii+spacing);
            previous2 = InterpolateValue(Test);
            sum += (previous+previous2);
            previous = previous2;
        }
        sum = sum*spacing*Point::Distance(p1,p2)/2;
        return sum;
        //    }
    }
    double Density::CalculateWeightedArea(const Poly Test) const{
        if (Values.empty()){
            throw std::runtime_error("Values have not been set!");
        }
        double sum = 0;
        int index = 0;
        double minx1,maxx1,miny1,maxy1,y0 = miny,x0 = minx;
        int size = Test.GetNVertices(), sizex = Ny-1;
        std::vector<bool> GridInPolyTemp = GetGridInRegion();
        if (size == 0){
            return Volume_Lower_Bound;
        }else{
            Test.GetExtrema(minx1, miny1, maxx1, maxy1);
            for(int ii=0;ii<Nx;ii++){
                if(x0<minx1 || x0>maxx1){
                    for(int jj=0;jj<Ny;jj++){
                        GridInPolyTemp[ii*Ny+jj] = false;
                    }
                }else{
                    y0 = miny;
                    for(int jj=0;jj<Ny;jj++){
                        if(y0<miny1 || y0>maxy1){
                            GridInPolyTemp[ii*Ny+jj] = false;
                        }else if(Test.pnpoly(Point(x0, y0))){
                            GridInPolyTemp[ii*Ny+jj] = true;
                            if (ii>0 && jj>0){
                                if (GridInPolyTemp[ii*Ny+jj] &&  GridInPolyTemp[(ii-1)*Ny+jj] && GridInPolyTemp[(ii-1)*Ny+jj-1] && GridInPolyTemp[ii*Ny+jj-1]){
                                    index = sizex*(ii-1)+jj-1;
                                    sum += Integral.Int[index];
                                }
                            }
                        }else{
                            GridInPolyTemp[ii*Ny+jj] = false;
                        }
                        y0+=dy;
                    }
                }
                x0+=dx;
            }
            
            if (sum>=Volume_Lower_Bound){
                return sum;
            }else{
                return Volume_Lower_Bound;
            }
        }
    }
    Point Density::CalculateCentroid(const Poly Test, const double &Volume) const{
        if (Values.empty()){
            throw std::runtime_error("Values have not been set!");
        }
        double sumx = 0, sumy = 0;
        int index = 0;
        double minx1,maxx1,miny1,maxy1,y0 = miny,x0 = minx;
        int size = Test.GetNVertices(), sizex = Ny-1;
        std::vector<bool> GridInPolyTemp = GetGridInRegion();
        if (size == 0){
            return Point();
        }else{
            Test.GetExtrema(minx1, miny1, maxx1, maxy1);
            for(int ii=0;ii<Nx;ii++){
                if(x0<minx1 || x0>maxx1){
                    for(int jj=0;jj<Ny;jj++){
                        GridInPolyTemp[ii*Ny+jj] = false;
                    }
                }else{
                    y0 = miny;
                    for(int jj=0;jj<Ny;jj++){
                        if(y0<miny1 || y0>maxy1){
                            GridInPolyTemp[ii*Ny+jj] = false;
                        }else if(Test.pnpoly(Point(x0, y0))){
                            GridInPolyTemp[ii*Ny+jj] = true;
                            if (ii>0 && jj>0){
                                if (GridInPolyTemp[ii*Ny+jj] &&  GridInPolyTemp[(ii-1)*Ny+jj] && GridInPolyTemp[(ii-1)*Ny+jj-1] && GridInPolyTemp[ii*Ny+jj-1]){
                                    index = sizex*(ii-1)+jj-1;
                                    sumx += Integral.Intx[index];
                                    sumy += Integral.Inty[index];
                                }
                            }
                        }else{
                            GridInPolyTemp[ii*Ny+jj] = false;
                        }
                        y0+=dy;
                    }
                }
                x0+=dx;
            }
            if (Volume<=Volume_Lower_Bound){
                return Point(minx1,miny1);
            }else{
                return Point(sumx/Volume, sumy/Volume);
            }
        }
    }
    // Partition Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    Partition::Partition(int NRegions, Density Prior, std::vector<double> desired_area, Parameters Alg_Params):NRegions(NRegions), Prior(Prior), desired_area(desired_area), Alg_Params(Alg_Params){Point::Robustness_Constant = Alg_Params.Robustness_Constant; CheckParams();}
    void Partition::SetPartitionVariables(int NRegions, Density Prior, std::vector<double> desired_area){this->NRegions = NRegions;this->Prior = Prior;this->desired_area = desired_area;CheckParams();}
    
    void Partition::CheckParams(){
        Prior.SetVolumeLowerBound(Alg_Params.Volume_Lower_Bound);
        double sum = 0;
        if (NRegions!=0 && desired_area.empty()){
            if (Alg_Params.Volume_Lower_Bound<=1.0/NRegions){
                std::vector<double> temp(NRegions, 1.0/NRegions);
                desired_area = temp;
            }else{
                throw std::runtime_error("Volume_Lower_Bound is too large for the number of regions. Try making Volume_Lower_Bound smaller or decreasing the number of regions");
            }
        }else if (desired_area.size() != NRegions){
            throw std::runtime_error("The size of desired_area must equal NRegions");
        }else{
            for (int ii = 0;ii<NRegions;ii++){
                if (desired_area[ii]<=Alg_Params.Volume_Lower_Bound){
                    throw std::runtime_error("Entries of desired_area must be greater than Alg_Params.Volume_Lower_Bound");
                }else{
                    sum+=desired_area[ii];
                }
            }
            if (sum!=1){
                std::cout<<"Warning: desired_areas not normalized; will be normalized automatically"<<std::endl;
                for (int ii = 0; ii<NRegions; ii++){
                    desired_area[ii]/=sum;
                    if (desired_area[ii]<Alg_Params.Volume_Lower_Bound){
                        throw std::runtime_error("Normalized areas too small. Decrease the number of regions, increase desired areas, or decrease Alg_Params.Volume_Lower_Bound");
                    }
                }
            }
            
        }
        
    }
    bool Partition::CreateDefaultCenters(const Poly Region, const double multiplier){
        Centers.clear();
        std::vector<Point> Vertices = Region.GetVertices();
        Point p1(Vertices[0]), p2(Vertices[1]), perp = Point::FindPerpDirection(p1, p2, multiplier);
        double spacing = 1.0/(NRegions+1);
        if (!Region.pnpoly(Point::FindPointAlongLine(p1, p2, 0.5).AddPoint(perp))){
            perp.FlipDirection();
        }
        for (int ii = 0; ii<NRegions;ii++){
            Centers.push_back(Point::FindPointAlongLine(p1, p2, spacing*(ii+1)).AddPoint(perp));
            if (!Region.pnpoly(Centers.back())){
                return false;
            }
        }
        return true;
        
    }
    void Partition::CreateDefaultCenters(const Poly Region, const double initial_multiplier, const int max_steps){
        double multiplier = initial_multiplier;
        for (int ii = 0;ii<max_steps+1;ii++){
            if (CreateDefaultCenters(Region, multiplier)){
                break;
            }else{
                if(ii == max_steps){
                    throw std::runtime_error("Unable to Create Default Centers");
                }else{
                    multiplier/=2;
                }
            }
        }
    }
    void Partition::InitializePartition(std::vector<Point> Centers, std::vector<double> Weights){
        Poly Region = Prior.GetRegion();
        std::vector<Poly> temp(NRegions);
        
        if (NRegions!=0 && Centers.empty()){
            int max_steps = 10;
            double initial_multiplier = 10e-3;
            CreateDefaultCenters(Region, initial_multiplier, max_steps);
        }else if(Centers.size() != NRegions){
            throw std::runtime_error("Centers must be the same size as NRegions!");
        }else{
            for(std::vector<Point>::iterator it = Centers.begin(); it!=Centers.end(); it++){
                if (!Region.pnpoly(*it)){
                    throw std::runtime_error("Centers must be located inside the region of interest");
                }
            }
            this->Centers = Centers;
        }
        
        if (NRegions!=0 && Weights.empty()){
            
            for (int ii = 0; ii<NRegions;ii++){
                this->Weights.push_back(0);
            }
        }else if (Weights.size()!= NRegions){
            throw std::runtime_error("Weights must be the same size as NRegions");
        }else{
            this->Weights = Weights;
        }
        Covering = temp;
    }
    void Partition::CreatePowerDiagram(void){
        double minx, maxx,miny, maxy = 0;
        Prior.GetExtrema(minx, miny, maxx, maxy);
        Poly Region = Prior.GetRegion();
        int NPoly = Region.GetNVertices(), count = 0;
        std::vector<Point> Vertices = Region.GetVertices();
        Point Test, p1, p2, p3, p4, p5, p6;
        double value1,previousvalue1 = 0, value2,previousvalue2 = 0, error;
        double tolerance = Alg_Params.Robustness_Constant, testvalue;
        long double increment;
        bool flag2 = false, whichpoly = false;
        long int mult = (int) 1/Alg_Params.Robustness_Constant;
        std::vector <Point> temp;
        Poly temp1;
        ClipperLib::Paths subj, clip(1), solution, temp2(1);
        ClipperLib::Clipper c;
        //Convert the base region into a format sutible for clipping
        for (int ii = 0;ii<NPoly;++ii){
            temp2[0].push_back(ClipperLib::IntPoint((long int)(Vertices[ii].x*mult),(long int)(Vertices[ii].y*mult)));
        }
        subj.push_back(temp2[0]);
        
        for (int ii = 0; ii<NRegions;++ii){
            solution = subj;
            for (int jj = 0; jj<NRegions;++jj){
                if (jj == ii){
                }else{
                    count = 0;
                    clip[0].clear();
                    value1 = 0;
                    testvalue = 0.50;
                    increment = 1;
                    flag2 = false;
                    error = INFINITY;
                    while (count<10000){
                        Test = Point::FindPointAlongLine(Centers[ii], Centers[jj], testvalue);
                        value1 = pow(Point::Distance(Centers[ii], Test),2) - Weights[ii];
                        value2 = pow(Point::Distance(Centers[jj], Test),2) - Weights[jj];
                        if (!flag2){
                            error = std::abs(value2-value1);
                            flag2 = true;
                            if (error<tolerance){
                                break;
                            }
                            previousvalue1 =value1;
                            previousvalue2 = value2;
                            if (value1>value2){
                                testvalue -=increment;
                            }else{
                                testvalue += increment;
                            }
                        }else{
                            error = std::abs(value2-value1);
                            if (error<tolerance){
                                break;
                            }
                            if ((value2>value1)&& (previousvalue2> previousvalue1)){
                                testvalue +=increment;
                                previousvalue1 = value1;
                                previousvalue2 = value2;
                            }else if ((value1>value2)&& (previousvalue1>previousvalue2)){
                                testvalue -=increment;
                                previousvalue1 = value1;
                                previousvalue2 = value2;
                            }else if (value2>value1 && previousvalue1>previousvalue2){
                                increment /= 2;
                                testvalue += increment;
                            }else {
                                increment /= 2;
                                testvalue -= increment;
                            }
                            count++;
                            
                        }
                    }
                    //Split the polygon into two sections
                    increment = 1;
                    while (true){
                        p1 = Test.AddPoint(Point::FindPerpDirection(Test, Centers[ii], Point::Distance(Test, Centers[ii])*increment));
                        p2 =Test.AddPoint(Point::FindPerpDirection(Test, Centers[ii], Point::Distance(Test, Centers[ii])*-increment));
                        if ((p1.x<minx && p2.x>maxx)||(p1.x>maxx && p2.x<minx)){
                            p3.x = p1.x;
                            p4.x = p2.x;
                            p5.x = p2.x;
                            p6.x = p1.x;
                            p3.y = std::min(miny,p1.y);
                            p3.y = std::min(p3.y,p2.y)-1;
                            p5.y = std::max(maxy,p1.y);
                            p5.y = std::max(p5.y,p2.y)+1;
                            p4.y = p3.y;
                            p6.y = p5.y;
                            break;
                        }else if ((p1.y<miny && p2.y>maxy)||(p1.y>maxy && p2.y<miny)){
                            p3.y = p1.y;
                            p4.y = p2.y;
                            p5.y = p2.y;
                            p6.y = p1.y;
                            p3.x = std::min(minx,p1.x);
                            p3.x = std::min(p3.x,p2.x)-1;
                            p5.x = std::max(maxx,p1.x);
                            p5.x = std::max(p5.x,p2.x)+1;
                            p4.x = p3.x;
                            p6.x = p5.x;
                            break;
                        }else {
                            increment = increment*2;
                        }
                        if (count >10000){
                            increment = increment * 500;
                        }
                    }
                    //Test to see which polygon is of relevance
                    temp = {p3, p4, p2, p1};
                    temp1.SetVertices(temp, false);
                    whichpoly = temp1.pnpoly(Centers[ii]);
                    if (whichpoly && (-Weights[ii]>pow(Point::Distance(Centers[ii],Centers[jj]),2)-Weights[jj])){
                        whichpoly = false;
                    }else if(!whichpoly && (-Weights[ii]>pow(Point::Distance(Centers[ii],Centers[jj]),2)-Weights[jj])){
                        whichpoly = true;
                    }
                    //Construct a polygon for clipping
                    if (whichpoly){
                        clip[0]<<ClipperLib::IntPoint((int)(p3.x*mult),(int)(p3.y*mult)) << ClipperLib::IntPoint((long int)(p4.x*mult),(long int)(p4.y*mult))<<ClipperLib::IntPoint((long int)(p2.x*mult),(long int)(p2.y*mult))<<ClipperLib::IntPoint((long int)(p1.x*mult),(long int)(p1.y*mult));
                        if (!ClipperLib::Orientation(clip[0])){
                            ClipperLib::ReversePath(clip[0]);
                        }
                    }else {
                        clip[0]<<ClipperLib::IntPoint((long int)(p5.x*mult),(long int)(p5.y*mult)) << ClipperLib::IntPoint((long int)(p6.x*mult),(long int)(p6.y*mult))<<ClipperLib::IntPoint((long int)(p1.x*mult),(long int)(p1.y*mult))<<ClipperLib::IntPoint((long int)(p2.x*mult),(long int)(p2.y*mult));
                        if (!ClipperLib::Orientation(clip[0])){
                            ClipperLib::ReversePath(clip[0]);
                        }
                    }
                    
                    //Clip the polygons
                    c.Clear();
                    c.AddPaths(solution,ClipperLib::ptSubject, true);
                    c.AddPaths(clip,ClipperLib::ptClip, true);
                    c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
                    for (int kk = 0;kk<solution.size();kk++){
                        ClipperLib::CleanPolygon(solution[0],1);
                    }
                }
                
            }
            //Convert solution to a Poly structure
            temp.clear();
            for (int pp = 0; pp<solution[0].size();++pp){
                temp.push_back(Point((double) solution[0][pp].X/(double)mult,(double) solution[0][pp].Y/(double)mult));
            }
            Covering[ii] = Poly(temp);
            
        }
        CleanCovering((double) 1.0/mult, mult);
        
    }
    
    void Partition::CleanCovering(const double tolerance, const long int &mult){
        double distance = 0;
        Poly temp;
        ClipperLib::Paths c(NRegions);
        std::vector<Point> Vert_ii, Vert_jj;
        for (int ii = 0;ii<NRegions;ii++){
            Vert_ii = Covering[ii].GetVertices();
            for (int jj = ii;jj<NRegions;jj++){
                Vert_jj = Covering[jj].GetVertices();
                for (int kk = 0;kk<Covering[ii].GetNVertices();kk++){
                    for (int pp = 0;pp<Covering[jj].GetNVertices();pp++){
                        distance = Point::Distance(Vert_ii[kk],Vert_jj[pp]);
                        if (distance<tolerance){
                            Vert_jj[pp] = Vert_ii[kk];
                        }
                    }
                }
                Covering[jj].SetVertices(Vert_jj, false);
            }
        }
        for (int ii = 0;ii<NRegions;ii++){
            Vert_ii = Covering[ii].GetVertices();
            for (int kk = 0;kk<Covering[ii].GetNVertices();kk++){
                
                c[ii]<<ClipperLib::IntPoint((long int)(Vert_ii[kk].x*mult),(long int)(Vert_ii[kk].y*mult));
            }
            
        }
        ClipperLib::CleanPolygons(c);
        for (int ii = 0;ii<NRegions;ii++){
            Vert_ii.clear();
            for (int kk = 0;kk<c[ii].size();++kk){
                Vert_ii.push_back(Point((double) c[ii][kk].X/(double)mult,(double) c[ii][kk].Y/(double)mult));
            }
            Covering[ii].SetVertices(Vert_ii);
        }
    }
    void Partition::CreateDelaunayGraph(DelaunayGraph &Delaunay) const{
        if (Delaunay.NRegions!=NRegions){
            throw std::runtime_error("DelaunayGraph has inconsistent sizes");
        }
        //DelaunayGraph Delaunay(NRegions);
        int NPoly_ii, NPoly_jj;
        std::vector<Point> Vert_ii, Vert_jj, temp;
        Point pi1, pi2, pj1, pj2, p;
        
        //Initialize Multi-Dimensional Array;
        bool flag = false;
        for (int ii = 0; ii<NRegions;ii++){
            NPoly_ii = Covering[ii].GetNVertices();
            Vert_ii = Covering[ii].GetVertices();
            for (int jj = ii+1;jj<NRegions;jj++){
                NPoly_jj = Covering[jj].GetNVertices();
                Vert_jj = Covering[jj].GetVertices();
                flag = false;
                for (int kk = 0;kk<NPoly_ii;kk++){
                    pi1 = Point(Vert_ii[kk].x,Vert_ii[kk].y);
                    if(flag){
                        break;
                    }
                    if (kk<NPoly_ii-1){
                        pi2 = Point(Vert_ii[kk+1].x,Vert_ii[kk+1].y);
                    }else{
                        pi2 = Point(Vert_ii[0].x,Vert_ii[0].y);
                    }
                    for (int pp = 0;pp<NPoly_jj;pp++){
                        pj1 = Point(Vert_jj[pp].x,Vert_jj[pp].y);
                        if (pp<NPoly_jj-1){
                            pj2 = Point(Vert_jj[pp+1].x,Vert_jj[pp+1].y);
                        }else{
                            pj2 = Point(Vert_jj[0].x,Vert_jj[0].y);
                        }
                        
                        if (Point::AreCollinear(pi1, pi2, pj1) && Point::AreCollinear(pi1, pi2, pj2)){
                            temp = Point::FindCollinearIntersection(pi1, pi2, pj1, pj2);
                            if (temp.empty()){
                                Delaunay.Graph[ii][jj][0] = p;
                                Delaunay.Graph[ii][jj][1] = p;
                                Delaunay.Graph[jj][ii][0] = Delaunay.Graph[ii][jj][0];
                                Delaunay.Graph[jj][ii][1] = Delaunay.Graph[ii][jj][1];
                            }else if (temp.size() == 1){
                                Delaunay.Graph[ii][jj][0] = temp[0];
                                Delaunay.Graph[ii][jj][1] = p;
                                Delaunay.Graph[jj][ii][0] = Delaunay.Graph[ii][jj][0];
                                Delaunay.Graph[jj][ii][1] = Delaunay.Graph[ii][jj][1];
                                
                            }else{
                                Delaunay.Graph[ii][jj][0] = temp[0];
                                Delaunay.Graph[ii][jj][1] = temp[1];
                                Delaunay.Graph[jj][ii][0] = Delaunay.Graph[ii][jj][0];
                                Delaunay.Graph[jj][ii][1] = Delaunay.Graph[ii][jj][1];
                                
                                flag = true;
                                break;
                            }
                            
                        }
                    }
                }
            }
        }
    }
    double Partition::GradientStepCenter(const std::vector<double> &volumes){
        Point Center, Center_ii, Errorxy;
        double Error = 0;
        for (int ii = 0; ii<NRegions;ii++){
            Center_ii = Centers[ii];
            Center_ii.FlipDirection();
            Center = Prior.CalculateCentroid(Covering[ii], volumes[ii]);
            Errorxy = Point::AddPoints(Center, Center_ii);
            Error+=Errorxy.Norm();
            Errorxy.Mult(Alg_Params.centers_step);
            Centers[ii] = Centers[ii].AddPoint(Errorxy);
            
        }
        
        return Error;
    }
    void Partition::GradientStepCenter(const double &temp_step, const std::vector<double> &volumes){
        if (temp_step <=0 || temp_step>1){
            throw std::runtime_error("temp_step must be between 0 and 1 (can be equal to 1, but not zero)");
        }
        Point Center;
        for (int ii = 0; ii<NRegions;ii++){
            Center = Prior.CalculateCentroid(Covering[ii], volumes[ii]);
            Centers[ii]=Point::FindPointAlongLine(Centers[ii], Center, temp_step);
        }
    }
    void Partition::GradientStepWeights(const std::vector<double> &volumes, const DelaunayGraph &SharedEdges){
        std::vector<double> totals(NRegions);
        Mult_Array integrals(NRegions), dist(NRegions), InvAreaDist(NRegions);
        
        for (int ii = 0; ii<NRegions;ii++){
            totals[ii] = 0;
            integrals.Array[ii][ii] = 0;
            dist.Array[ii][ii] = 0;
            InvAreaDist.Array[ii][ii] = 0;
            for (int jj = ii+1;jj<NRegions;jj++){
                dist.Array[ii][jj] = Point::Distance(Centers[ii], Centers[jj]);
                dist.Array[jj][ii] = dist.Array[ii][jj];
                InvAreaDist.Array[ii][jj] = (desired_area[jj]/volumes[jj])-(desired_area[ii]/volumes[ii]);
                InvAreaDist.Array[jj][ii] = -InvAreaDist.Array[ii][jj];
                if (!isinf(SharedEdges.Graph[ii][jj][1].x)){
                    integrals.Array[ii][jj] = Prior.LineIntegral(Alg_Params.line_int_step, SharedEdges.Graph[ii][jj][0], SharedEdges.Graph[ii][jj][1]);
                }else{
                    integrals.Array[ii][jj] = 0;
                }
                integrals.Array[jj][ii] = integrals.Array[ii][jj];
                totals[ii] = totals[ii]+InvAreaDist.Array[ii][jj]*(1/dist.Array[ii][jj])*integrals.Array[ii][jj];
            }
            for (int jj=0;jj<ii;jj++){
                totals[ii] = totals[ii]+InvAreaDist.Array[ii][jj]*(1/dist.Array[ii][jj])*integrals.Array[ii][jj];
                
            }
            
            if(Covering[ii].GetNVertices() == 0){
                Weights[ii] +=Alg_Params.weights_step*2;
            }else{
                Weights[ii] += - totals[ii]*Alg_Params.weights_step;
            }
        }
        
    }
    double Partition::CalculateError(const std::vector<double> &volumes){
        double sum = 0;
        for (int jj = 0; jj<NRegions;jj++){
            sum+=std::abs(volumes[jj]-desired_area[jj])*std::abs(volumes[jj]-desired_area[jj]);
        }
        return sum;
    }
    
    std::vector<double> Partition::CalculateVolumes(void){
        std::vector<double> result(NRegions);
        for (int jj = 0;jj<NRegions;jj++){
            result[jj] = Prior.CalculateWeightedArea(Covering[jj]);
        }
        return result;
    }
    void Partition::CalculatePartition(bool WriteToFile, std::string filename_partition, std::string filename_centers){
        if (Prior.GetRegion().GetNVertices() == 0){
            throw std::runtime_error("Prior has not been initialized");
        }else if (Centers.empty()){
            throw std::runtime_error("Centers and Weights have not been initialized");
        }
        std::ofstream file1, file2;
        int count1, count2;
        std::vector<double> volumes(NRegions);
        std::vector<Point> Vert;
        DelaunayGraph Delaunay(NRegions);
        double initial_step = 1, error = INFINITY, error_vol = INFINITY;
        if (WriteToFile){
            file1.open(filename_centers);
            file2.open(filename_partition);
        }
        if (WriteToFile){
            for (int ii = 0; ii<NRegions;ii++){
                file1<<Centers[ii].x<<","<<Centers[ii].y<<std::endl;
                Vert = Covering[ii].GetVertices();
                for (int jj = 0; jj<Covering[ii].GetNVertices();jj++){
                    file2<<Vert[jj].x<<","<<Vert[jj].y<<" ";
                }
                file2<<std::endl;
            }
            file1<<std::endl;
            file2<<std::endl;
        }
        CreatePowerDiagram();
        if (WriteToFile){
            for (int ii = 0; ii<NRegions;ii++){
                file1<<Centers[ii].x<<","<<Centers[ii].y<<std::endl;
                Vert = Covering[ii].GetVertices();
                for (int jj = 0; jj<Covering[ii].GetNVertices();jj++){
                    file2<<Vert[jj].x<<","<<Vert[jj].y<<" ";
                }
                file2<<std::endl;
            }
            file1<<std::endl;
            file2<<std::endl;
        }
        volumes = CalculateVolumes();
        GradientStepCenter(initial_step, volumes);
        CreatePowerDiagram();
        if (WriteToFile){
            for (int ii = 0; ii<NRegions;ii++){
                file1<<Centers[ii].x<<","<<Centers[ii].y<<std::endl;
                Vert = Covering[ii].GetVertices();
                for (int jj = 0; jj<Covering[ii].GetNVertices();jj++){
                    file2<<Vert[jj].x<<","<<Vert[jj].y<<" ";
                }
                file2<<std::endl;
            }
            file1<<std::endl;
            file2<<std::endl;
        }
        count2 = 0;
        while (error>Alg_Params.convergence_criterion && count2<Alg_Params.max_iterations_centers){
            volumes = CalculateVolumes();
            error_vol = CalculateError(volumes);
            
            std::cout<<error_vol<<std::endl;
            count1 = 0;
            while (error_vol >Alg_Params.volume_tolerance &&count1<Alg_Params.max_iterations_volume){
                CreateDelaunayGraph(Delaunay);
                GradientStepWeights(volumes, Delaunay);
                CreatePowerDiagram();
                if (WriteToFile){
                    for (int ii = 0; ii<NRegions;ii++){
                        file1<<Centers[ii].x<<","<<Centers[ii].y<<std::endl;
                        Vert = Covering[ii].GetVertices();
                        for (int jj = 0; jj<Covering[ii].GetNVertices();jj++){
                            file2<<Vert[jj].x<<","<<Vert[jj].y<<" ";
                        }
                        file2<<std::endl;
                    }
                    file1<<std::endl;
                    file2<<std::endl;
                }
                volumes = CalculateVolumes();
                error_vol = CalculateError(volumes);
                std::cout<<count1<<std::endl;
                std::cout<<error_vol<<std::endl;
                count1++;
            }
            error = GradientStepCenter(volumes);
            std::cout<<error<<std::endl;
            CreatePowerDiagram();
            if (WriteToFile){
                for (int ii = 0; ii<NRegions;ii++){
                    file1<<Centers[ii].x<<","<<Centers[ii].y<<std::endl;
                    Vert = Covering[ii].GetVertices();
                    for (int jj = 0; jj<Covering[ii].GetNVertices();jj++){
                        file2<<Vert[jj].x<<","<<Vert[jj].y<<" ";
                    }
                    file2<<std::endl;
                }
                file1<<std::endl;
                file2<<std::endl;
            }
            count2++;
        }
        if (WriteToFile){
            for (int ii = 0; ii<NRegions;ii++){
                file1<<Centers[ii].x<<","<<Centers[ii].y<<std::endl;
                Vert = Covering[ii].GetVertices();
                for (int jj = 0; jj<Covering[ii].GetNVertices();jj++){
                    file2<<Vert[jj].x<<","<<Vert[jj].y<<" ";
                }
                file2<<std::endl;
            }
            file1<<std::endl;
            file2<<std::endl;
        }
        
    } 
}
