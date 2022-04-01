#include <cmath>
#include <math.h>
#include <iostream>

double func(double x){
	return (sin(x)*pow(x,3));
	//return pow(x,sin(x))+5;
}

double fib(int n){
	double first = pow((1.0+sqrt(5))/2,n);
	double second = pow((1.0-sqrt(5))/2,n);
	return (first-second)/sqrt(5);
}

double dichotomy(double left, double right, double eps){
	double sigma = eps/2-eps/1000;
	double x1, x2, f1, f2;
	int n = 0, count = 0;
	while(fabs(right - left) > eps){
		n++;
		x1 = (left+right)/2-sigma;
		x2 = (left+right)/2+sigma;
		f1 = func(x1);
		f2 = func(x2);
		count+=2;
		if(f1 < f2)
			right=x2;
		else
			left=x1;
	}
	std::cout<<"Dichotomy method:"<<std::endl;
	std::cout<<"Iterations: "<<n<<" Function calls: "<<count<<std::endl;
	std::cout<<"Result: "<<(left+right)/2<<" "<<"Y = "<<func((left+right)/2)<<std::endl;
	return (left+right)/2;
}

double goldenSection(double left, double right, double eps){
    double resphi = (3-sqrt(5))/2;
    double x1 = left+resphi*(right-left);
    double x2 = right-resphi*(right-left);
    double f1 = func(x1);
    double f2 = func(x2);
    int n = 0, count = 2;
    while (fabs(right - left) > eps){
    	n++;
    	if (f1<f2)
		{
	        right = x2;
	        x2 = x1;
	        f2 = f1;
	        x1 = left+resphi*(right-left);
	        f1 = func(x1);
	        count++;
		}
      	else
		{
	        left = x1;
	        x1 = x2;
	        f1 = f2;
	        x2 = right-resphi*(right-left);
	        f2 = func(x2);
	        count++;
		}
	}
	std::cout<<"Golden Section method:"<<std::endl;
	std::cout<<"Iterations: "<<n<<" Function calls: "<<count-1<<std::endl;
	std::cout<<"Result: "<<(x1+x2)/2<<" "<<"Y = "<<func((x1+x2)/2)<<std::endl;
    return (x1+x2)/2;
}

double fibonacci(double left, double right, double eps){
	int n = 1;
	while((right-left)/eps>fib(n+2))
		n++;
	double x1 = left+fib(n)*(right-left)/fib(n+2);
	double x2 = left+fib(n+1)*(right-left)/fib(n+2);
	double f1 = func(x1);
    double f2 = func(x2);
	int count = 2;
	for(int k=1;k<=n;k++){
    	if (f1<f2)
		{
	        right = x2;
	        x2 = x1;
	        f2 = f1;
	        x1 = left+fib(n-k+1)*(right-left)/fib(n-k+3);
	        f1 = func(x1);
	        count++;
		}
      	else
		{
	        left = x1;
	        x1 = x2;
	        f1 = f2;
	        x2 = left+fib(n-k+2)*(right-left)/fib(n-k+3);
	        f2 = func(x2);
	        count++;
		}
	}
	std::cout<<"Fibonacci method:"<<std::endl;
	std::cout<<"Iterations: "<<n<<" Function calls: "<<count-1<<std::endl;
	std::cout<<"Result: "<<(x1+x2)/2<<" "<<"Y = "<<func((x1+x2)/2)<<std::endl;
    return (x1+x2)/2;
}

double parabolas(double left, double right, double eps){
	double x = (left+right)/2;
	double f1 = func(left);
	double f2 = func(x);
	double f3 = func(right);
	double u = x-(pow((x-left),2)*(x-left)*(f2-f3)-pow((x-right),2)*(f2-f1))/(2*((x-left)*(f2-f3)-(x-right)*(f2-f1)));
	double fu = func(u);
	int n = 0, count = 4;
	while(fabs(right - left) > eps){
		n++;
		if(f2>fu){
			f3 = f2;
			right = x;
			x = u;
			f2 = fu;
			u = x-(pow((x-left),2)*(f2-f3)-pow((x-right),2)*(f2-f1))/(2*((x-left)*(f2-f3)-(x-right)*(f2-f1)));
			fu=func(u);
			count++;
		}
		else{
			f1 = f2;
			left = x;
			x = u;
			f2 = fu;
			u = x-(pow((x-left),2)*(f2-f3)-pow((x-right),2)*(f2-f1))/(2*((x-left)*(f2-f3)-(x-right)*(f2-f1)));
			fu=func(u);
			count++;
		}
	}
	std::cout<<"Parabolas method:"<<std::endl;
	std::cout<<"Iterations: "<<n<<" Function calls: "<<count-1<<std::endl;
	std::cout<<"Result: "<<(left+right)/2<<" "<<"Y = "<<func((left+right)/2)<<std::endl;
    return (left+right)/2;
}

double combinedBrent(double left, double right, double eps){
	double k = (3-sqrt(5))/2;
	double x = (left+right)/2, om = x, v = x;
	double fx = func(x), fom = fx, fv = fx;
	double d = right-left, e = d, g, u, fu;
	int n = 0, count = 1;
	while(fabs(right - left) > eps){
		if(n>=1000)
			break;
		n++;
		g = e;
		e = d;
		if(x!=om&&x!=v&&v!=om&&fx!=fom&&fx!=fv&&fv!=fom){
			if(om>v)
				u = x-(pow((x-v),2)*(x-v)*(fx-fom)-pow((x-om),2)*(fx-fv))/(2*((x-v)*(fx-fom)-(x-om)*(fx-fv)));
			else
				u = x-(pow((x-om),2)*(x-om)*(fx-fv)-pow((x-v),2)*(fx-fom))/(2*((x-om)*(fx-fv)-(x-v)*(fx-fom)));
		}
		if(u>=(left+eps)&&u<=(right-eps)&&fabs(u-x)<g/2)
			d = fabs(u-x);
		else{
			if(x<(right-left)/2){
				u = x+k*(right-x);
				d = right-x;
			}
			else{
				u = x-k*(x-left);
				d = x-left;
			}
		}
		if(fabs(u-x)<eps){
			if(u-x>=0)
				u = x+eps;
			else
				u = x-eps;
		}
		fu=func(u);
		count++;
		if(fu<=fx){
			if(u>=x){
				left = x;
			}
			else{
				right = x;
			}
			v = om;
			om = x;
			x = u;
			fv = fom;
			fom = fx;
			fx = fu;
		}
		else{
			if(u>=x){
				right = u;
			}
			else{
				left = x;
			}
			if(fu<=fom||om==x){
				v = om;
				om = u;
				fv = fom;
				fom = fu;
			}
			else if(fu<=fv||v==x||v==om){
				v=u;
				fv=fu;
			}
		}
	}
	std::cout<<"Combined Brent method:"<<std::endl;
	std::cout<<"Iterations: "<<n<<" Function calls: "<<count-1<<std::endl;
	std::cout<<"Result: "<<(left+right)/2<<" "<<"Y = "<<func((left+right)/2)<<std::endl;
    return (left+right)/2;
}

int main(){
	double left = -3, right = 3;
	for(int i=1; i<6; i++){
		double eps = pow(0.1,i);
		std::cout<<std::endl;
		std::cout<<"Current epsilon = "<<eps<<std::endl;
		std::cout<<std::endl<<std::endl;
		dichotomy(left, right, eps);
		std::cout<<std::endl;
		goldenSection(left, right, eps);
		std::cout<<std::endl;
		fibonacci(left, right, eps);
		std::cout<<std::endl;
		parabolas(left, right, eps);
		std::cout<<std::endl;
		combinedBrent(left, right, eps);
		std::cout<<std::endl;
	}
}
