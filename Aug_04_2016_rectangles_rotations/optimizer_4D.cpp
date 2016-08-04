// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    optimization routines from the dlib C++ Library.

    The library provides implementations of the conjugate gradient,  BFGS,
    L-BFGS, and BOBYQA optimization algorithms.  These algorithms allow you to
    find the minimum of a function of many input variables.  This example walks
    though a few of the ways you might put these routines to use.


		This is a useful video to setup dlib in visual studio:

		https://www.youtube.com/watch?v=lhfIbaygA-s

		and a website related with the video:

		http://quantlabs.net/blog/2014/04/how-to-get-great-open-source-dlib-c-bayesian-network-library-working-on-visual-studio/#sthash.9OJNS8p6.dpuf

*/
// OpenCV includes
#include <opencv2/opencv.hpp>



#include <dlib/optimization.h>
#include <iostream>


#include "functions.h"
#include "Utilities.h"
#include <vector>

#include "simple_optimizer.h"
#include <iostream>




using namespace std;
//using namespace dlib;




// ----------------------------------------------------------------------------------------

// In dlib, the general purpose solvers optimize functions that take a column
// vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.

//typedef matrix<double,0,1> column_vector;

typedef matrix<double,0,1> column_vector;


// Parameters of our sigmoid-like 
//double mu = 0.5;
// Class for functions
functions getFunct;

//double xCoord = double(35);
//double yCoord = double(35);
//double r = double(15);
//double mu = 0.5;
 
//int xMin = 30;
//int yMin = 30;
//int rMin = 5;
//double muMin = 0.1;


//int xMax = 70;
//int yMax = 70;
//int rMax = 20;
//double muMax = 1.0;



// ----------------------------------------------------------------------------------------
// Below we create a few functions.  When you get down into main() you will see that
// we can use the optimization algorithms to find the minimums of these functions.
// ----------------------------------------------------------------------------------------

double rosen (const column_vector& m)
/*
    This function computes what is known as Rosenbrock's function.  It is 
    a function of two input variables and has a global minimum at (1,1).
    So when we use this function to test out the optimization algorithms
    we will see that the minimum found is indeed at the point (1,1). 
*/
{

		// arguments of the function
    const double x = m(0); 
    const double y = m(1);
		const double r = m(2);
		const double mu = m(3);


		//return -1.0*getFunct.ObjFunction4D(x,y,r,mu,lambdaVal,rPower);
		return 0;
}


// This is a helper function used while optimizing the rosen() function.  
const column_vector rosen_derivative (const column_vector& m)
/*!
    ensures
        - returns the gradient vector for the rosen function
!*/
{
    const double x = m(0);
    const double y = m(1);
		const double r = m(2);
		const double mu = m(3);

		Point4D gradVal;


    // make us a column vector of length 3
    column_vector res(4);


		//gradVal = getFunct.ObjFunction4DGrad(x,y,r,mu,lambdaVal,rPower);
		
		res(0) = -1.0*gradVal[0];
		res(1) = -1.0*gradVal[1];
		res(2) = -1.0*gradVal[2];
		res(3) = -1.0*gradVal[3];

    
		return res;
}




// This function computes the Hessian matrix for the rosen() fuction.  This is
// the matrix of second derivatives.
matrix<double> rosen_hessian (const column_vector& m)
{ 
    const double x = m(0);
    const double y = m(1);
		const double r = m(2);
		const double mu = m(3);

    matrix<double> res(4,4);


    // now compute the second derivatives 
    
		//res(0,0) = 1200*x*x - 400*y + 2; // second derivative with respect to x
    //res(1,0) = res(0,1) = -400*x;   // derivative with respect to x and y
    //res(1,1) = 200;                 // second derivative with respect to y
    

		matrix<double> hessVal(4,4);

		//hessVal = getFunct.ObjFunction3DHessian(x,y,r,mu);

		//res = -1.0*getFunct.ObjFunction4DHessian(x,y,r,mu,lambdaVal,rPower);

		///////////////////////////////////////////////////////
		//
		//
		//             SILLY TEST TO CHECK HESSIAN
		//
		//
		//////////////////////////////////////////////////////


		//res(0,0) = -1.0;
		//res(1,1) = -1.0;
		//res(2,2) = -1.0;
		//res(3,3) = -1.0;
		//res(0,1) = res(1,0) = 0.0;
		//res(0,2) = res(0,2) = 0.0;
		//res(0,3) = res(0,3) = 0.0;
		//res(1,2) = res(2,1) = 0.0;
		//res(1,3) = res(3,1) = 0.0;
		//res(2,3) = res(3,2) = 0.0;



		//for (int i=0; i<=3 ; i++)
		//{
		//	res(i,i) = hessVal(i,i) ;
		//	for (int j=i+1; j<=3; j++)
		//	{
		//		res(i,j) = res(j,i) = hessVal(i,j);
		//	}
		//}

		return res;
}

// ----------------------------------------------------------------------------------------

class test_function
{
    /*
        This object is an example of what is known as a "function object" in C++.
        It is simply an object with an overloaded operator().  This means it can 
        be used in a way that is similar to a normal C function.  The interesting
        thing about this sort of function is that it can have state.  
        
        In this example, our test_function object contains a column_vector 
        as its state and it computes the mean squared error between this 
        stored column_vector and the arguments to its operator() function.

        This is a very simple function, however, in general you could compute
        any function you wanted here.  An example of a typical use would be 
        to find the parameters of some regression function that minimized 
        the mean squared error on a set of data.  In this case the arguments
        to the operator() function would be the parameters of your regression
        function.  You would loop over all your data samples and compute the output 
        of the regression function for each data sample given the parameters and 
        return a measure of the total error.   The dlib optimization functions 
        could then be used to find the parameters that minimized the error.
    */
public:

    test_function (
        const column_vector& input
    )
    {
        target = input;
    }

    double operator() ( const column_vector& arg) const
    {
        // return the mean squared error between the target vector and the input vector
        return mean(squared(target-arg));
    }

private:
    column_vector target;
};

// ----------------------------------------------------------------------------------------

class rosen_model 
{
    /*!
        This object is a "function model" which can be used with the
        find_min_trust_region() routine.  
    !*/

public:
    typedef ::column_vector column_vector;
    typedef matrix<double> general_matrix;

    double operator() (
        const column_vector& x
    ) const { return rosen(x); }

    void get_derivative_and_hessian (
        const column_vector& x,
        column_vector& der,
        general_matrix& hess
    ) const
    {
        der = rosen_derivative(x);
        hess = rosen_hessian(x);
    }
};

// ----------------------------------------------------------------------------------------

int main()
{
	

	//cout << "Let's test dlib idiot!";
	//cout << endl;
	//cin.get();

	//imshow("Culero", Test);
	//cv::waitKey();
	//exit(0);


	//for (int i=0;i<1000;i++)
	//	pTest[i+50000] = 1.0;


	//for (int j = 250000; j<750000;j++)
	//	for (int i = 250; i < 750; i++)
	//		pTest[i + j] = 1.0;




	//cv::circle(Test,cv::Point(500,500),100,0.5,-1);
	//cv::rectangle(Test,cv::Rect(200,200,100,100),0.5,1);
	//cv::imshow("puto",Test);
	//cv::waitKey();
	//exit(0);

	simple_optimizer myGrad;

	// Data parameters
	int xInt = 200;
	int yInt = 200;  
	int radiusInt = 50;
	double muInit = 0.8;
	double theta = 27.15;


	// Synthetic image
	cv::Mat Test(1000,1000,CV_64F,0.0);
	double * pTest = &Test.at<double>(0,0);
	
	// RECTANGLE
	cv::rectangle(Test,cv::Point(xInt-radiusInt/0.8,yInt - radiusInt/0.8),
		cv::Point(xInt+radiusInt/0.8,yInt + radiusInt/0.8),1.0,-1);
	
	// Reference for rotations:
	// http://felix.abecassis.me/2011/10/opencv-rotation-deskewing/
	// http://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine#warpaffine
	cv::Mat rot_mat =  getRotationMatrix2D(cv::Point(xInt,yInt),theta,1.0);
	//cv::Mat rotated;
	
	// Apply rotation...
	cv::warpAffine(Test, Test, rot_mat, Test.size());
	//cv::imshow("Rotated?" , Test);
	//cv::waitKey(10);
	//cin.get();




	// CIRCLE
	//cv::circle(Test,cv::Point(xInt,yInt),radiusInt,0.5,-1);


	// ROTATED RECTANGLE

	//cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(xInt,yInt),
	//	cv::Size2f(radiusInt,2.0*radiusInt),25.0);
	//cv::Point2f vertices[4];
	//rRect.points(vertices);
	//for (int i = 0; i < 4; i++){
		//line(Temp, vertices[i], vertices[(i+1)%4], cv::Scalar(0,255,0));
	//	line(Test, vertices[i], vertices[(i+1)%4], 0.5);
	//}



	// Set the synthetic image as the input for the gradient descent fn.
	myGrad.SetImageData(Test);


	// Initial guess
	double xCoord;
	double yCoord;
	double r;
	double mu;

	// for regularization...
	int rPower = 2;
	double lambdaVal = 0.0; //2.0;

	xCoord = 350;
	yCoord = 350;
	r = 80;
	//mu = 0.5;
	mu = 0.05; // 0.05;


	cout << " Initial coordinates ..." << endl;
	cout << xInt << " " << " " << yInt << " " << radiusInt << " " << muInit << endl;

	myGrad.gradDescent(xCoord, yCoord, r, mu, lambdaVal,rPower);


	cin.get();

	//exit(0);




	// Rectangle

	int uppXint = 100;
	int uppYint = 100;
	int lowXint = 200;
	int lowYint = 200;



	// For constrained search...


	int xMin = 60;
	int yMin = 60;
	int rMin = 40;
	double muMin = 0.004;


	int xMax = 1000;
	int yMax = 1000;
	int rMax = 500;
	double muMax = 2.0;


	
	//exit(0);
	
	
	
	
	
	//char test_Alg = 'all';
  // options for algorithm
	// test_Alg = 'all' -- run all algorithms ---> DEFAULT
	//					= 'bfgs' -- bfgs with analytical gradient
	//					= 'bfgs_num' -- bfgs with numerical gradient
	//					= 'lbfgs' -- lbfgs with analytical gradient
	//					= 'lbfgs_num' -- lbfgs with numerical gradient
	//					= 'lbfgs_box' -- lbfgs with analytical gradient constrained to a box
	//					= 'newton' -- Newton search analytical Hessian 
	//					= 'newton_box' -- Newton search analytical Hessian constrained to a box
	//          = 'trust_radius' -- Trust radius




    try
    {



        // make a column vector of length 3
        column_vector starting_point(4);

				// set limits of box
				column_vector boxMin(4), boxMax(4);

				boxMin(0) = double(xMin);
				boxMin(1) = double(yMin);
				boxMin(2) = double(rMin);
				boxMin(3) = muMin;


				boxMax(0) = double(xMax);
				boxMax(1) = double(yMax);
				boxMax(2) = double(rMax);
				boxMax(3) = muMax;


        // Set the starting point to (4,8).  This is the point the optimization algorithm
        // will start out from and it will move it closer and closer to the function's 
        // minimum point.   So generally you want to try and compute a good guess that is
        // somewhat near the actual optimum value.
        //starting_point = 28, 28;

				starting_point = xCoord, yCoord, r, mu;





        // The first example below finds the minimum of the rosen() function and uses the
        // analytical derivative computed by rosen_derivative().  Since it is very easy to
        // make a mistake while coding a function like rosen_derivative() it is a good idea
        // to compare your derivative function against a numerical approximation and see if
        // the results are similar.  If they are very different then you probably made a 
        // mistake.  So the first thing we do is compare the results at a test point:

        //cout << "Difference between analytic derivative and numerical approximation of derivative: " 
        //      << length(derivative(rosen)(starting_point) - rosen_derivative(starting_point)) << endl;



				cout << "Computing derivatives: " 
					<< endl;
				cout << " Numerical  " << endl;
				cout << derivative(rosen)(starting_point) << endl;
				cout << " Analytical " << endl;  
				cout << rosen_derivative(starting_point) << endl;

				cout << "Difference between analytic derivative and numerical approximation of derivative: " 
					<< endl;
				cout << length(derivative(rosen)(starting_point) - rosen_derivative(starting_point)) << endl;
				cout << endl;


				//cin.get();




				//cv::Mat Test2(1000,1000,CV_8UC3,0.0);
				//cv::circle(Test2,cv::Point(500,500),100,cv::Scalar(255,0,0),-1);
				//cv::rectangle(Test2,cv::Rect(200,200,100,100),cv::Scalar(255,255,255),1);
				//cv::imshow("puto2",Test2);
				//cv::waitKey();


			
				//cout << "***********************************************************" << endl;
        //cout << "Find the minimum of the function" << endl;
				//cout << "Using BFGS algorithm..." << endl;
				//cout << "It uses analytic gradient of function ... " <<endl;

				//cin.get();


        // Now we use the find_min() function to find the minimum point.  The first argument
        // to this routine is the search strategy we want to use.  The second argument is the 
        // stopping strategy.  Below I'm using the objective_delta_stop_strategy which just 
        // says that the search should stop when the change in the function being optimized 
        // is small enough.

        // The other arguments to find_min() are the function to be minimized, its derivative, 
        // then the starting point, and the last is an acceptable minimum value of the rosen() 
        // function.  That is, if the algorithm finds any inputs to rosen() that gives an output 
        // value <= -1 then it will stop immediately.  Usually you supply a number smaller than 
        // the actual global minimum.  So since the smallest output of the rosen function is 0 
        // we just put -1 here which effectively causes this last argument to be disregarded.

        //find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
        //         objective_delta_stop_strategy(1e-7), // Stop when the change in rosen() is less than 1e-7
        //         rosen, rosen_derivative, starting_point, -3000);
        // Once the function ends the starting_point vector will contain the optimum point 
        // of (1,1).
        //cout << "rosen solution:\n" << starting_point << endl;



				//cin.get();




				//cout << "***********************************************************" << endl;
				//cout << "Find the minimum of the function" << endl;
				//cout << "Using BFGS algorithm..." << endl;
				//cout << "It uses approximate gradient (numerical) of function ... " <<endl;


        // Now let's try doing it again with a different starting point and the version
        // of find_min() that doesn't require you to supply a derivative function.  
        // This version will compute a numerical approximation of the derivative since 
        // we didn't supply one to it.
        //starting_point = -94, 5.2;

				//starting_point = xCoord, yCoord, r, mu;


        //find_min_using_approximate_derivatives(bfgs_search_strategy(),
        //                                       objective_delta_stop_strategy(1e-7),
        //                                       rosen, starting_point, -3000);
        // Again the correct minimum point is found and stored in starting_point
        //cout << "rosen solution:\n" << starting_point << endl;

				

				//cin.get();



				// abort , exit code
				//abort();
				//exit(EXIT_FAILURE);




				//cout << "***********************************************************" << endl;
				//cout << "Find the minimum of the function" << endl;
				//cout << "Using L-BFGS algorithm..." << endl;
				//cout << "It uses analytic gradient of function ... " <<endl;


        // Here we repeat the same thing as above but this time using the L-BFGS 
        // algorithm.  L-BFGS is very similar to the BFGS algorithm, however, BFGS 
        // uses O(N^2) memory where N is the size of the starting_point vector.  
        // The L-BFGS algorithm however uses only O(N) memory.  So if you have a 
        // function of a huge number of variables the L-BFGS algorithm is probably 
        // a better choice.
        //starting_point = 0.8, 1.3;
        
				//starting_point = xCoord, yCoord, r, mu;

				//find_min(lbfgs_search_strategy(10),  // The 10 here is basically a measure of how much memory L-BFGS will use.
        //         objective_delta_stop_strategy(1e-7),  // Adding be_verbose() causes a message to be 
                                                                    // printed for each iteration of optimization.
       //          rosen, rosen_derivative, starting_point, -3000);

        //cout << endl << "rosen solution: \n" << starting_point << endl;


				//cin.get();



        //starting_point = -94, 5.2;




				//cout << "***********************************************************" << endl;
				//cout << "Find the minimum of the function" << endl;
				//cout << "Using L-BFGS algorithm..." << endl;
				//cout << "It uses approximate gradient (numerical) of function ... " <<endl;

        

				//starting_point = xCoord, yCoord, r, mu;
				//find_min_using_approximate_derivatives(lbfgs_search_strategy(10),
        //                                       objective_delta_stop_strategy(1e-7),
        //                                       rosen, starting_point, -3000);
        //cout << "rosen solution: \n"<< starting_point << endl;



				//cin.get();



				// VERY IMPORTANT!!!
				// OPTIMIZATION WITH BOUNDS CONTRAINTS -- WITHIN A BOX
				// SEEMS TO BE WHAT WE NEED!!!


        // dlib also supports solving functions subject to bounds constraints on
        // the variables.  So for example, if you wanted to find the minimizer
        // of the rosen function where both input variables were in the range
        // 0.1 to 0.8 you would do it like this:


        //starting_point = 0.1, 0.1; // Start with a valid point inside the constraint box.
        //find_min_box_constrained(lbfgs_search_strategy(10),  
        //                         objective_delta_stop_strategy(1e-9),  
        //                         rosen, rosen_derivative, starting_point, 0.1, 0.8);



				// Initial drawing...
				cv::Mat Test(1000,1000,CV_8UC3,0.0);
				cv::rectangle(Test,cv::Rect(uppXint,uppYint,lowXint,lowYint),cv::Scalar(255,255,255),-1);
				cv::circle(Test,cv::Point(xInt,yInt),radiusInt,cv::Scalar(255,0,0),mu);
				cv::imshow("Initial gradient",Test);
				cv::waitKey();


				// To print elements of the matrix...
				//cout << "M = "<< endl << " "  << Test << endl << endl;




				cout << "***********************************************************" << endl;
				cout << "Find the minimum of the function CONSTRAINED TO A BOX" << endl;
				cout << "Using L-BFGS algorithm..." << endl;
				cout << "It uses analytical gradient of function ... " <<endl;


				


				starting_point = xCoord, yCoord, r, mu; // Start with a valid point inside the constraint box.
				cout << find_min_box_constrained(lbfgs_search_strategy(10),
																 objective_delta_stop_strategy(1e-7),  
																	rosen, rosen_derivative, starting_point, boxMin, boxMax);


				cout << endl;

				cin.get();


        // Here we put the same [0.1 0.8] range constraint on each variable, however, you
        // can put different bounds on each variable by passing in column vectors of
        // constraints for the last two arguments rather than scalars.  

        //cout << endl << "constrained rosen solution: \n" << starting_point << endl;

				cout << " Final coordinates ..." << endl;
				cout << starting_point(0) << " " << starting_point(1) << " " << \
					starting_point(2) << " " << starting_point(3) <<  endl;

				//xInt = int( starting_point(0) );
				//yInt = int( starting_point(1) );
				//radiusInt = int( starting_point(2) );
				//mu = starting_point(3);

				// Final drawing...
				cv::Mat Test2(1000,1000,CV_8UC3,0.0);
				cv::rectangle(Test2,cv::Rect(uppXint,uppYint,lowXint,lowYint),cv::Scalar(255,255,255),-1);
				cv::circle(Test2,cv::Point(int( starting_point(0) ),int( starting_point(1) )),
					int( starting_point(2) ),cv::Scalar(255,0,0),starting_point(3));
				cv::imshow("Final gradient",Test2);
				cv::waitKey();
				

				
				cin.get();
				exit(0);

				//exit(EXIT_FAILURE);





				//////////////////////////////////////////////////////////////
				//
				//
				//             CHECK THE HESSIAN YOU IDIOT !!!
				//
				//
				/////////////////////////////////////////////////////////////

        // In many cases, it is useful if we also provide second derivative information
        // to the optimizers.  Two examples of how we can do that are shown below.  
        //starting_point = 0.8, 1.3;
        




				// Initial drawing...
				cv::Mat Test3(1000,1000,CV_8UC3,0.0);
				cv::rectangle(Test3,cv::Rect(uppXint,uppYint,lowXint,lowYint),cv::Scalar(255,255,255),-1);
				cv::circle(Test3,cv::Point(xInt,yInt),radiusInt,cv::Scalar(255,0,0),mu);
				cv::imshow("Initial Hessian",Test3);
				cv::waitKey();




				cout << "***********************************************************" << endl;
				cout << "Find the minimum of the function..." << endl;
				cout << "Using Newton search algorithm..." << endl;
				cout << "It uses analytical Hessian of function ... " <<endl;


				starting_point = double(xInt),double(yInt), double(radiusInt), muInit;
				
				//cout << "With Hessian? " <<  endl;
				
				find_min(newton_search_strategy(rosen_hessian),
                 //objective_delta_stop_strategy(1e-7),
								 objective_delta_stop_strategy(1e-7),
                 rosen,
                 rosen_derivative,
                 starting_point,
                 -3000.0);
        cout << "rosen solution: \n"<< starting_point << endl;


				cin.get();





				cout << " Final coordinates ..." << endl;
				cout << starting_point(0) << " " << starting_point(1) << " " << \
					starting_point(2) << " " << starting_point(3) <<  endl;

				xInt = int( starting_point(0) );
				yInt = int( starting_point(1) );
				radiusInt = int( starting_point(2) );
				mu = starting_point(3);

				// Final drawing...
				cv::Mat Test4(1000,1000,CV_8UC3,0.0);
				cv::rectangle(Test4,cv::Rect(uppXint,uppYint,lowXint,lowYint),cv::Scalar(255,255,255),-1);
				cv::circle(Test4,cv::Point(int( starting_point(0) ),int( starting_point(1) )),
					int( starting_point(2) ),cv::Scalar(255,0,0),starting_point(3));
				cv::imshow("Final Hessian",Test4);
				cv::waitKey();

				cin.get();

				//////////////////////////////////////////////////////////////
				//
				//
				//             CHECK THE HESSIAN YOU IDIOT !!!
				//
				//
				/////////////////////////////////////////////////////////////

        // We can also use find_min_trust_region(), which is also a method which uses
        // second derivatives.  For some kinds of non-convex function it may be more
        // reliable than using a newton_search_strategy with find_min().
        //starting_point = 0.8, 1.3;



				//cout << "***********************************************************" << endl;
				//cout << "Find the minimum of the function CONSTRAINED TO A BOX" << endl;
				//cout << "Using Newton search..." << endl;
				//cout << "It uses analytical Hessian of function ... " <<endl;



				//////////////////////////////////////////////////////////////
				//
				//
				//             CHECK THE HESSIAN YOU IDIOT !!!
				//						 Note the convergence criterion has been relaxed!
				//
				//
				/////////////////////////////////////////////////////////////



				//starting_point = xCoord, yCoord, r, mu; // Start with a valid point inside the constraint box.
				//cout << find_min_box_constrained(newton_search_strategy(rosen_hessian),
				//	objective_delta_stop_strategy(1e-2),  
				//	rosen, rosen_derivative, starting_point, boxMin, boxMax);




				//////////////////////////////////////////////////////////////
				//
				//
				//             CHECK THE HESSIAN YOU IDIOT !!!
				//
				//
				/////////////////////////////////////////////////////////////



				//cout << endl;


				// Here we put the same [0.1 0.8] range constraint on each variable, however, you
				// can put different bounds on each variable by passing in column vectors of
				// constraints for the last two arguments rather than scalars.  

				//cout << endl << "constrained rosen solution with Hessian: \n" << starting_point << endl;

				//cin.get();


				exit(EXIT_FAILURE);




				cout << "***********************************************************" << endl;
				cout << "Find the minimum of the function with trust radius" << endl;

				//						 Noticed the convergence criterion has been relaxed!
				
				starting_point = xCoord, yCoord, r, mu;


        find_min_trust_region(objective_delta_stop_strategy(1e-5), //.be_verbose(),
            rosen_model(), 
            starting_point, 
            2 // initial trust region radius
        );
        cout << "rosen solution: \n"<< starting_point << endl;

				

				cout << "Done!" << endl;


				cin.get();







				//////////////////////////////////////////////////////////////


				//                  STOP CODE HERE!!!


				//////////////////////////////////////////////////////////////


				//exit(EXIT_FAILURE);


    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;


				cin.get();
    }

		


}

