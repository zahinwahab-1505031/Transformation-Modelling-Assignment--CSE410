#include <bits/stdc++.h>

using namespace std;
#define PI acos(-1.0)
struct Point
{
    double x;
    double y;
    double z;
    double w;
    Point()
    {
        x = 0;
        y = 0;
        z = 0;
        w = 1;
    }
    Point(double x1,double y1,double z1)
    {
        x = x1;
        y = y1;
        z = z1;
        w = 1;
    }
    void normalize()
    {
        double d = sqrt(x*x+y*y+z*z);
        x/=d;
        y/=d;
        z/=d;
    }
    void printPoint()
    {
        cout << "("<< x << "," << y << "," << z << ")" << endl;
    }
    Point multiply(double a)
    {
        Point ret;
        ret.x = x*a;
        ret.y = y*a;
        ret.z = z*a;
        return ret;
    }
};
Point crossProduct(Point p,Point q)
{
    /*************
    |x  y   z|
    |x1 y1 z1|
    |x2 y2 z2|

    x(y1z2-z1y2)
    y(z1x2 - x1z2)
    z(x1y2 - x2y1)



    ******************/
    Point cp;
    cp.x = p.y*q.z - p.z*q.y;
    cp.y = p.z*q.x - p.x*q.z;
    cp.z = p.x*q.y - p.y*q.x;
//    cout << "After cross product ";
//    p.printPoint();
//    q.printPoint();
//    cout << "we get: \n";
//    cp.printPoint();
    return cp;
}
Point add(Point p,Point q)
{
    /*************
    ******************/
    Point cp;
    cp.x = p.x+q.x;
    cp.y = p.y+q.y;
    cp.z = p.z+q.z;
    cp.w =1;

//    cp.x /= cp.w;
//    cp.y /= cp.w;
//    cp.z /= cp.w;
//    cp.w /= cp.w;
//    cout << "After adding ";
//    p.printPoint();
//    q.printPoint();
//    cout << "we get: \n";
//    cp.printPoint();

    return cp;
}
Point subtract(Point p,Point q)
{
    /*************
    ******************/
    Point cp;
    cp.x = p.x-q.x;
    cp.y = p.y-q.y;
    cp.z = p.z-q.z;
    cp.w = 0;

//    cp.x /= cp.w;
//    cp.y /= cp.w;
//    cp.z /= cp.w;
//    cp.w /= cp.w;
//    cout << "After adding ";
//    p.printPoint();
//    q.printPoint();
//    cout << "we get: \n";
//    cp.printPoint();

    return cp;
}
double dotProduct(Point p,Point q)
{
    /*************
    |x  y   z|
    |x1 y1 z1|
    |x2 y2 z2|

    x(y1z2-z1y2)
    y(z1x2 - x1z2)
    z(x1y2 - x2y1)



    ******************/
    double cp = p.x*q.x + p.y*q.y+p.z*q.z;
    return cp;
}
struct Matrix_41
{
    double element[4];
    Matrix_41()
    {
        for(int i=0; i<4; i++) element[i] = 0;
    }
    Matrix_41(double given[4])
    {
        for(int i=0; i<4; i++) element[i] = given[i];
    }
    void printColumnVector()
    {
        cout << "("<< element[0] << ","<< element[1] << "," << element[2] << "," << element[3] << ")" << endl;
    }

};
struct Matrix
{
    Matrix_41 rowVector[4];
    Matrix()
    {
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                if(i==j)  rowVector[i].element[j] = 1;
                else rowVector[i].element[j] = 0;
            }
        }
    }
    void printMatrix()
    {
        cout << "--------------------\n";
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                cout << rowVector[i].element[j]<< " ";
            }
            cout << endl;
        }
        cout << "---------------------\n";
    }
    void scale(double a)
    {
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                rowVector[i].element[j]/=a;
            }

        }
    }
};
Matrix matrixMultiply(Matrix matrix_1, Matrix matrix_2)
{
    Matrix ret;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double sum = 0;
            for (int k = 0; k < 4; k++)
            {
                sum+= matrix_1.rowVector[i].element[k] * matrix_2.rowVector[k].element[j];

            }
            ret.rowVector[i].element[j] = sum;
        }
    }




    return ret;
}
void PrintStack(stack<int> s)
{
    if(s.empty()) return;
    int x = s.top();

    s.pop();
    PrintStack(s);

    cout << x << " ";
    s.push(x);
}
void PrintStack(stack<pair<Matrix,int> > s)
{
    if(s.empty()) return;
    pair <Matrix,int> p = s.top();

    s.pop();
    PrintStack(s);

    cout << p.second << " ";
    s.push(p);
}
Point eye,look,up;
double fovY,aspectRatio,near,far;
stack <pair<Matrix,int> > Stack;
stack <int> instructionStack;
double degreeToRadian(double angle){
    return (PI*angle/180.0);
}
void projectionTransformation()
{
    ifstream inFile;
    inFile.open("stage2.txt");
    ofstream outFile;
    outFile.open("stage3.txt");
    if (!inFile)
    {
        cerr << "Unable to open file stage2.txt";
        exit(1);
    }

    Matrix V;
    double fovX = fovY * aspectRatio;

    double t = near * tan(degreeToRadian(fovY/2));
    double r = near * tan(degreeToRadian(fovX/2));
    /*************
    near/r 0 0 0
     0 near/t 0 0
      0 0 -(far+near)/(far-near) -(2*far*near)/(far-near)
       0 0 -1 0
    *********/
    for(int i=0;i<4;i++) V.rowVector[i].element[i] = 0;
    V.rowVector[0].element[0] = near/r;
    V.rowVector[1].element[1] = near/t;
    V.rowVector[2].element[2] = -(far+near)/(far-near);
    V.rowVector[2].element[3]= -(2*far*near)/(far-near);
    V.rowVector[3].element[2] = -1;
    //V.printMatrix();



    int cnt=0;
        double var1;
    while(inFile>>var1)
    {
        cnt++;
        double Y[4][1];
        Y[0][0] = var1;


        inFile >> Y[1][0] >> Y[2][0];
        Y[3][0] = 1;

        double Z[4][1];
        for (int i =0 ; i< 4; i++ )      //first matrix row
        {
            for (int j = 0; j < 1; j++ )          //second matrix column
            {
                double sum = 0;
                for (int k = 0; k < 4; k++ )      //second matrix row==first matrix column
                {
                    sum += V.rowVector[i].element[k] * Y[k][j];
                }
                Z[i][j] = sum;
            }
        }
         for(int i=0; i<3; i++)
            for(int j=0; j<1; j++)
                    Z[i][j] /= Z[3][j];
        for(int i=0; i<3; i++)
            for(int j=0; j<1; j++)
                   // cout << Z[j][i] << " ";
               outFile << setprecision(7)<<fixed<< Z[j][i] << " ";
        outFile << endl;
       if(cnt%3==0) outFile << endl;
       // cout << endl;
    }



}
void viewTransformation()
{
    ifstream inFile;
    inFile.open("stage1.txt");
    ofstream outFile;
    outFile.open("stage2.txt");
    if (!inFile)
    {
        cerr << "Unable to open file stage1.txt";
        exit(1);
    }
    double p;
//    while(inFile>>p){
//
//    }
    Point l = subtract(look,eye);
    l.normalize();
    Point r = crossProduct(l,up);
    r.normalize();
    Point u = crossProduct(r,l);
//    cout << "l: ";
//    l.printPoint();
//    cout << "r: ";
//    r.printPoint();
//    cout << "u: ";
//    u.printPoint();
    Matrix T;
    T.rowVector[0].element[3] = -eye.x;
    T.rowVector[1].element[3] = -eye.y;
    T.rowVector[2].element[3] = -eye.z;
//    cout << "T:\n";
//    T.printMatrix();
    Matrix R;
    R.rowVector[0].element[0] = r.x;
    R.rowVector[0].element[1] = r.y;
    R.rowVector[0].element[2] = r.z;

    R.rowVector[1].element[0] = u.x;
    R.rowVector[1].element[1] = u.y;
    R.rowVector[1].element[2] = u.z;

    R.rowVector[2].element[0] = -l.x;
    R.rowVector[2].element[1] = -l.y;
    R.rowVector[2].element[2] = -l.z;
//    cout << "R:\n";
//    R.printMatrix();
    Matrix V = matrixMultiply(R,T);
//    cout << "V:\n";
//    V.printMatrix();
    int cnt=0;
    double var1;
    while(inFile>>var1)
    {
        cnt++;
        double Y[4][1];
        Y[0][0] = var1;

        inFile >> Y[1][0] >> Y[2][0];
        Y[3][0] = 1;

        double Z[4][1];
        for (int i =0 ; i< 4; i++ )      //first matrix row
        {
            for (int j = 0; j < 1; j++ )          //second matrix column
            {
                double sum = 0;
                for (int k = 0; k < 4; k++ )      //second matrix row==first matrix column
                {
                    sum += V.rowVector[i].element[k] * Y[k][j];
                }
                Z[i][j] = sum;
            }
        }

        for(int i=0; i<3; i++)
            for(int j=0; j<1; j++)
                outFile << setprecision(7)<<fixed<< Z[j][i] << " ";
        outFile << endl;
        if(cnt%3==0) outFile  << endl;
    }

    outFile.close();
    projectionTransformation();






}
int main()
{

    Matrix identityMatrix;
    //identityMatrix.printMatrix();
    Stack.push(make_pair(identityMatrix,0));
    ifstream inFile;
    inFile.open("scene.txt");
    ofstream outfile;
    outfile.open ("stage1.txt");
    if (!inFile)
    {
        cerr << "Unable to open file scene.txt";
        exit(1);
    }
    string p;
    int var = 0;
    int ins = 0;
    while(inFile>>p)
    {

        if(var==0)eye.x = atof(p.c_str());
        else if(var==1) eye.y = atof(p.c_str());
        else if(var==2) eye.z = atof(p.c_str());
        else if(var==3) look.x = atof(p.c_str());
        else if(var==4) look.y = atof(p.c_str());
        else if(var==5) look.z = atof(p.c_str());
        else if(var==6) up.x = atof(p.c_str());
        else if(var==7) up.y = atof(p.c_str());
        else if(var==8) up.z = atof(p.c_str());
        else if(var==9) fovY = atof(p.c_str());
        else if(var==10) aspectRatio = atof(p.c_str());
        else if(var==11) near = atof(p.c_str());
        else if(var==12) far = atof(p.c_str());
        var++;
        if(p=="triangle")
        {

            //cout << "TRIANGLE"<<endl;
            Matrix_41 columnVector[3];
            for(int i=0; i<3; i++)
            {
                for(int j=0; j<4; j++)
                {
                    if(j<3) inFile >> columnVector[i].element[j];
                    else  columnVector[i].element[j]=1;
                }
            }
//            for(int i=0; i<3; i++)
//            {
//                columnVector[i].printColumnVector();
//            }
            double Y[4][3];
            for(int i=0; i<4; i++)
            {
                for(int j=0; j<3; j++)
                {
                    Y[i][j] = columnVector[j].element[i];

                }

            }
//            cout << "Y VECTOR: \n";
//            for(int i=0; i<4; i++)
//            {
//                for(int j=0; j<3; j++)
//                {
//                    cout << Y[i][j] << " ";
//
//                }
//                cout << endl;
//
//            }

            Matrix top = Stack.top().first;
//            cout << "STACK HEAD IS \n";
//            top.printMatrix();
            Point result;
            double sum = 0;
            double Z[4][3];
            for (int i =0 ; i< 4; i++ )      //first matrix row
            {
                for (int j = 0; j < 3; j++ )          //second matrix column
                {
                    sum = 0;
                    for (int k = 0; k < 4; k++ )      //second matrix row==first matrix column
                    {
                        sum += top.rowVector[i].element[k] * Y[k][j];
                    }
                    Z[i][j] = sum;
                }
            }
//            cout << "Z VECTOR: \n";
//            for(int i=0; i<4; i++)
//            {
//                for(int j=0; j<3; j++)
//                {
//                    cout << Z[i][j] << " ";
//
//                }
//                cout << endl;
//
//            }
            for(int i=0; i<3; i++)
            {
                for(int j=0; j<3; j++)
                {
                    outfile << setprecision(7)<<fixed<< Z[j][i] << " ";

                }
                outfile << endl;

            }
            outfile << endl;




        }
        else if(p=="translate")
        {
            ins++;
//            cout << "========================INS: "<< ins <<"====================="<< endl;
//            cout << "TRANSLATE"<<endl;
            Matrix translationMatrix;
            inFile >> translationMatrix.rowVector[0].element[3] >> translationMatrix.rowVector[1].element[3] >>
                   translationMatrix.rowVector[2].element[3];
            // translationMatrix.printMatrix();
            Matrix top = Stack.top().first;
//            cout  << "TOP: "<< endl;
//            top.printMatrix();
//            cout  << "PUSHED MATRIX: "<< endl;

            Matrix res = matrixMultiply(top,translationMatrix);
            //res.scale(res.rowVector[3].element[3]);
            //res.printMatrix();
            Stack.push(make_pair(res,ins));

        }
        else if(p=="rotate")
        {
            ins++;

//            cout << "========================INS: "<< ins <<"====================="<< endl;
//            cout << "ROTATE" <<endl;
            Matrix rotationMatrix;
            Point a;
            double angle;
            inFile >> angle >> a.x >> a.y >> a.z;
            //a.printPoint();
            a.normalize();
            //a.printPoint();
            double cos_t = cos(PI*angle/180.0);
            double sin_t = sin(PI*angle/180.0);
            Point c1,c2,c3;
            Point i(1,0,0),j(0,1,0),k(0,0,1);
            Point res1,res2,res3;
            res1 = crossProduct(a,i);
            res1 = res1.multiply(sin_t);
            res2 = crossProduct(a,j);
            res2 = res2.multiply(sin_t);
            res3 = crossProduct(a,k);
            res3 = res3.multiply(sin_t);
            res1 = add(i.multiply(cos_t),res1);
            res2 = add(j.multiply(cos_t),res2);
            res3 = add(k.multiply(cos_t),res3);
            Point res11,res22,res33;
            res11 = a.multiply(dotProduct(a,i));
            res11 = res11.multiply(1-cos_t);
            res22 = a.multiply(dotProduct(a,j));
            res22 = res22.multiply(1-cos_t);
            res33 = a.multiply(dotProduct(a,k));
            res33 = res33.multiply(1-cos_t);

            c1 = add(res1,res11);
            c2 = add(res2,res22);
            c3 = add(res3,res33);
//            cout << "C1: " << endl;
//            c1.printPoint();
//            cout << "C2: " << endl;
//            c2.printPoint();
//            cout << "C3: " << endl;
//            c3.printPoint();
//            rotate(angle,a.x,a.y,a.z);
            /***
            c1.x c2.x c3.x 0
            c1.y c2.y c3.y 0
            c1.z c2.z c3.z 0
            0 0 0 1 ***/
            rotationMatrix.rowVector[0].element[0] = c1.x;
            rotationMatrix.rowVector[0].element[1] = c2.x;
            rotationMatrix.rowVector[0].element[2] = c3.x;
            rotationMatrix.rowVector[1].element[0] = c1.y;
            rotationMatrix.rowVector[1].element[1] = c2.y;
            rotationMatrix.rowVector[1].element[2] = c3.y;
            rotationMatrix.rowVector[2].element[0] = c1.z;
            rotationMatrix.rowVector[2].element[1] = c2.z;
            rotationMatrix.rowVector[2].element[2] = c3.z;
            //rotationMatrix.printMatrix();
            Matrix top = Stack.top().first;
            Matrix res = matrixMultiply(top,rotationMatrix);
            //res.scale(res.rowVector[3].element[3]);
            Stack.push(make_pair(res,ins));

        }
        else if(p=="scale")
        {
            ins++;
//            cout << "========================INS: "<< ins <<"====================="<< endl;
            Matrix scalingMatrix;
            double sx,sy,sz;
            inFile >> sx >> sy >> sz;
            scalingMatrix.rowVector[0].element[0] = sx;
            scalingMatrix.rowVector[1].element[1] = sy;
            scalingMatrix.rowVector[2].element[2] = sz;
//            cout << "SCALE\n";
//            scalingMatrix.printMatrix();
            Matrix top = Stack.top().first;
//            cout << "Top matrix is:\n";
//            top.printMatrix();
            if(Stack.size()==1) Stack.push(make_pair(scalingMatrix,ins));
            else
            {
                Matrix res = matrixMultiply(top,scalingMatrix);
                //res.scale(res.rowVector[3].element[3]);
//                cout << "Multiplying with top matrix we get:\n";
//                res.printMatrix();
                Stack.push(make_pair(res,ins));
            }

        }
        else if(p=="push")
        {

            ins++;

//            cout << "========================INS: "<< ins <<"====================="<< endl;
//            cout << "PUSH\n";
            instructionStack.push(ins);
//            cout << "Printing instruction stack\n";
//            PrintStack(instructionStack);
//            cout << endl;
//            cout << "After pushing Printing tranformation stack\n";
//            PrintStack(Stack);
//            cout << endl;
        }
        else if(p=="pop")
        {

            ins++;

//            cout << "========================INS: "<< ins <<"====================="<< endl;
//            cout << "POP\n";
            int n = instructionStack.top();
            instructionStack.pop();
//            cout << "Before popping Printing tranformation stack\n";
//            PrintStack(Stack);
//            cout << endl;
//            cout << "Printing instruction stack\n";
//            PrintStack(instructionStack);
//            cout << endl;
            while(!Stack.empty())
            {
                pair <Matrix,int> var = Stack.top();
                int x = var.second;
                if(x>n)
                {
                    Stack.pop();


                }
                else break;
            }
//            cout << "After popping Printing tranformation stack\n";
//            PrintStack(Stack);
//            cout << endl;

        }



    }

//    eye.printPoint();
//    look.printPoint();
//    up.printPoint();
     outfile.close();

    viewTransformation();


    return 0;


}
