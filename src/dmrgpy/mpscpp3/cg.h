// conjugate gradient method to solve A*x = b for MPS

static auto solveAxb=[](auto A, auto b )
{
   double TOLERANCE = 1.0e-10;
   auto X = b*0.; // initialize
   auto R = b;
   auto P = R;
   int maxm = get_int_value("maxm") ; // bond dimension for KPM
   auto cutoff = get_float_value("cutoff") ; // cutoff for KPM
   int k = 0; // counter
   auto pp = {"MaxDim",maxm,"Cutoff",cutoff} ; // DMRG parameters
   auto n = 100 ; // maximum number of iterations
   while ( k < n )
   {
      auto Rold = R;   
      auto AP = applyMPO(A,P,pp) ;
      double alpha = inner( prime(R), R ) / max(inner(prime(P), AP ), NEARZERO );
      X = sum( X, alpha*P,pp );   
      R = sum( R, -alpha*AP,pp );  
      if ( inner(prime(R),R ) < TOLERANCE ) break; 
      double beta = inner(prime(R), R ) / max(inner(prime(Rold), Rold ), NEARZERO );
      P = vectorCombination( 1.0, R, beta, P );
      k++;
   }
   return X;
};




// conjugate gradient method to solve [(H-w)**2 +d**2]*x = b for MPS

static auto solveHwdxb=[](auto H, auto w, auto d, auto b )
{
   double TOLERANCE = 1.0e-10;
   auto X = b*0.; // initialize
   auto R = b;
   auto P = R;
   int maxm = get_int_value("maxm") ; // bond dimension for KPM
   auto cutoff = get_float_value("cutoff") ; // cutoff for KPM
   int k = 0; // counter
   auto pp = {"MaxDim",maxm,"Cutoff",cutoff} ; // DMRG parameters
   auto n = 100 ; // maximum number of iterations
   while ( k < n )
   {
      auto Rold = R;   
// compute A*x
      auto Hx = applyMPO(H,P,pp) ; // Apply H*x
      auto H2x = applyMPO(H,Hx,pp) ; // Apply H^2*x
      auto AP = sum(H2x,-2*w*Hx,dd) ; // Add first contribution
      AP = sum(AP,(w**2+d**2)*P) ; // add last contribution
      double alpha = inner(prime(R), R ) / max(inner(prime(P), AP ), NEARZERO );
      X = sum( X, alpha*P,pp );   
      R = sum( R, -alpha*AP,pp );  
      if ( inner(prime(R),R ) < TOLERANCE ) break; 
      double beta = inner(prime(R), R ) / max(inner(prime(Rold), Rold ), NEARZERO );
      P = vectorCombination( 1.0, R, beta, P );
      k++;
   }
   return X;
}
;
