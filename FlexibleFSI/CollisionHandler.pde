/**********************************************************************
      CollisionHandler class: ControlPoints and/or FlexibleSheets
      are gathered by this class to handle any collisions during
      the simulation. Detection and resolution are separated.

Example code:
      (To be filled)

**********************************************************************/
class CollisionHandler {
  //================= Attributes ====================//
  int Ncp; // number of control points to consider
  int Nsp; // number of springs to consider
  
  ArrayList<ControlPoint> LocalCPoints = new ArrayList<ControlPoint>();
  ArrayList<Spring> LocalSprings = new ArrayList<Spring>();
  
  ArrayList<ControlPoint> FastCPi = new ArrayList<ControlPoint>();
  ArrayList<ControlPoint> FastCPj = new ArrayList<ControlPoint>();
  FloatList FastT = new FloatList();
  
  ArrayList<ControlPoint> CP = new ArrayList<ControlPoint>();
  ArrayList<Spring> SP = new ArrayList<Spring>();
  FloatList RewTcp = new FloatList();
  FloatList RewTp1 = new FloatList();
  FloatList RewTp2 = new FloatList();
  
  
  //================= Constructor ====================//
  CollisionHandler() {
    Ncp = 0;
    Nsp = 0;
    // Just exploit the class for its methods
  }
  CollisionHandler( ControlPoint [] cpoints ) {
    Ncp = cpoints.length;
    Nsp = 0;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
  }
  CollisionHandler( ControlPoint cpoint ) {
    Ncp = 1;
    Nsp = 0;
    LocalCPoints.add(cpoint);
  }
  CollisionHandler( ControlPoint [] cpoints, Spring [] springs ) {
    Ncp = cpoints.length;
    Nsp = springs.length;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
    for (int i=0; i<Nsp; i++) LocalSprings.add(springs[i]);
  }
  CollisionHandler( ControlPoint [] cpoints, Spring spring ) {
    Ncp = cpoints.length;
    Nsp = 1;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
    LocalSprings.add(spring);
  }
  
  
  //================= Method ====================//
  // Detect Boundary Collisions
  void DetectBoundCollision() {
    for (int i=0; i<Ncp; i++) {
      ControlPoint cp = LocalCPoints.get(i);
      if (cp.position.x < cp.diameter/2) {
        this.ResolveBoundCollisions( "West", cp );
      }
      else if (cp.position.x > cp.myWindow.x.inE - cp.diameter/2) {
        this.ResolveBoundCollisions( "East", cp );
      }
      if (cp.position.y < cp.diameter/2) {
        this.ResolveBoundCollisions( "North", cp );
      }
      else if (cp.position.y > cp.myWindow.y.inE - cp.diameter/2) {
        this.ResolveBoundCollisions( "South", cp );
      }
    }
  }
  // Detect Boundary Collisions
  boolean DetectBoundCollision( ControlPoint cp ) {
    boolean Flag = false;
    if (cp.position.x < cp.diameter/2) {
      this.ResolveBoundCollisions( "West", cp );
      Flag = true;
    }
    else if (cp.position.x > cp.myWindow.x.inE - cp.diameter/2) {
      this.ResolveBoundCollisions( "East", cp );
      Flag = true;
    }
    if (cp.position.y < cp.diameter/2) {
      this.ResolveBoundCollisions( "North", cp );
      Flag = true;
    }
    else if (cp.position.y > cp.myWindow.y.inE - cp.diameter/2) {
      this.ResolveBoundCollisions( "South", cp );
      Flag = true;
    }
    return Flag;
  }
  // Resolve boundary collisions
  void ResolveBoundCollisions( String bound, ControlPoint cp ){
    if ((bound.equals("North")==true) || (bound.equals("South")==true)) {
      float vx = cp.velocity.x;
      float vy = -cp.velocity.y;
      float x = cp.position.x;
      float y = cp.position.y;
      if (bound.equals("North")==true) y = cp.diameter/2;
      else if (bound.equals("South")==true) y = cp.myWindow.y.inE-cp.diameter/2;
      cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
    }
    if ((bound.equals("West")==true) || (bound.equals("East")==true)) {
      float vx = -cp.velocity.x;
      float vy = cp.velocity.y;
      float x = cp.position.x;
      float y = cp.position.y;
      if (bound.equals("West")==true) x = cp.diameter/2;
      else if (bound.equals("East")==true) x = cp.myWindow.x.inE-cp.diameter/2;
      cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
    }
  }
  
  // Detect CPoint-CPoint collisions
  void DetectCPointCPointCollision() {
    int N = LocalCPoints.size();
    if (N>1) {
      for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
          ControlPoint cpi = LocalCPoints.get(i);
          ControlPoint cpj = LocalCPoints.get(j);
          float dij = cpi.distance(cpj);
          float clearRad = (cpi.diameter + cpj.diameter)/2.;
          if (dij<=clearRad) {
            float penet = 0.5*(cpi.diameter + cpj.diameter) - dij;
            float rewindi = penet/(0.5*cpi.diameter);
            float rewindj = penet/(0.5*cpj.diameter);
            this.ResolveCPCPCollisions( cpi, rewindi, cpj, rewindj );
            // Must include a detection method for fast moving cpoints
            // Cheking if their relative velocity is greater than clearRad 
          }
        }
      }
    }
  }
  // Detect CPoint-CPoint collisions
  boolean DetectCPointCPointCollision( ControlPoint cpi, ControlPoint cpj ) {
    boolean Flag = false;
    float dij = cpi.distance(cpj);
    float clearRad = (cpi.diameter + cpj.diameter)/2.;
    if (dij<=clearRad) {
      float penet = 0.5*(cpi.diameter + cpj.diameter) - dij;
      float rewindi = penet/(0.5*cpi.diameter);
      float rewindj = penet/(0.5*cpj.diameter);
      this.ResolveCPCPCollisions( cpi, rewindi, cpj, rewindj );
      Flag = true;
      // Must include a detection method for fast moving cpoints
     // Cheking if their relative velocity is greater than clearRad 
    }
    return Flag;
  }
  // Resolve CPoint-CPoint collisions
  void ResolveCPCPCollisions( ControlPoint cpi,  float rewti, ControlPoint cpj, float rewtj ) {
    cpi.rewindPosition(rewti);
    cpj.rewindPosition(rewtj);
    float Vxi = (cpi.velocity.x*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.x;
    float Vyi = (cpi.velocity.y*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.y;
    float Vxj = (cpj.velocity.x*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.x;
    float Vyj = (cpj.velocity.y*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.y;
    cpi.UpdateVelocity( Vxi, Vyi );
    cpj.UpdateVelocity( Vxj, Vyj );
  }
  
  // Detect Fast moving CPoint-CPoint collisions
  void DetectFastCPCPCollision() {
    int N = LocalCPoints.size();
    if (N>1) {
      for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
          ControlPoint cpi = LocalCPoints.get(i);
          ControlPoint cpj = LocalCPoints.get(j);
          float tol = 1e-7;
          float tcol, denom, Rt;
          float vxi, vyi, vxj, vyj;
          float xoldi, yoldi, xoldj, yoldj;
          vxi = cpi.position.x-cpi.positionOld.x;
          vyi = cpi.position.y-cpi.positionOld.y;
          vxj = cpj.position.x-cpj.positionOld.x;
          vyj = cpj.position.y-cpj.positionOld.y;
          xoldi = cpi.positionOld.x; yoldi = cpi.positionOld.y;
          xoldj = cpj.positionOld.x; yoldj = cpj.positionOld.y;
          denom = sq(vxi-vxj)+sq(vyi-vyj);
          if (denom<tol) {
            // Moving in parallel or not moving at all....
          }
          else {
            tcol = -1*((vxi-vxj)*(xoldi-xoldj)+(vyi-vyj)*(yoldi-yoldj))/denom;
            Rt = sq((xoldi+tcol*vxi) - (xoldj + tcol*vxj)) + sq((yoldi+tcol*vyi) - (yoldj + tcol*vyj));
            float clearRad = (cpi.diameter + cpj.diameter)/2.;
            if (Rt<sq(clearRad)) {
              if ((tcol>=0) && (tcol<=1)) {
                // Fast collisions happened!
                FastCPi.add(cpi);
                FastCPj.add(cpj);
                FastT.append(tcol);
              }
            }
          }
        }
      }
    }
  }
  
  // Resolve Fast moving CPoint-CPoint collisions
  void ResolveFastCPCPCollisions() {
    int N = FastCPi.size();
    for (int ij=0; ij<N; ij++) {
      ControlPoint cpi = FastCPi.get(ij);
      ControlPoint cpj = FastCPj.get(ij);
      float tcol = FastT.get(ij);
      
      float Vxi = (cpi.velocity.x*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.x;
      float Vyi = (cpi.velocity.y*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.y;
      float Vxj = (cpj.velocity.x*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.x;
      float Vyj = (cpj.velocity.y*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.y;
      cpi.UpdateVelocity( Vxi, Vyi );
      cpj.UpdateVelocity( Vxj, Vyj );
      //cpi.impUpdate(tcol);
      //cpj.impUpdate(tcol);
      delay(3000);
    }
    FastCPi = new ArrayList<ControlPoint>();
    FastCPj = new ArrayList<ControlPoint>();
    FastT = new FloatList();
  }
  
  // Maybe check for parallel spring...?
  // Still some false-positive checks occuring
  // Detect point-spring collisions
  void  DetectCPointSpringCollision() {
    if ((Nsp>0) && (Ncp>2)) {
      for (int icp=0; icp<Ncp; icp++) {
        ControlPoint cp = LocalCPoints.get(icp);
        for (int isp=0; isp<Nsp; isp++) {
          Spring sp = LocalSprings.get(isp);
          ControlPoint p1, p2;
          p1 = sp.p1;
          p2 = sp.p2;
          if ((cp!=p1) && (cp!=p2)) {
            LineSweepsPoint( sp, cp );
          }
        }
      }
    }
  }
  // Detect point-spring collisions
  boolean  DetectCPointSpringCollision( Spring sp, ControlPoint cp ) {
    boolean Flag = false;
    ControlPoint p1, p2;
    p1 = sp.p1;
    p2 = sp.p2;
    if ((cp!=p1) && (cp!=p2)) {
      Flag = LineSweepsPoint( sp, cp );
    }
    return Flag;
  }
  // Check if a line sweeps a point
  boolean LineSweepsPoint( Spring sp, ControlPoint cp ) {
    boolean Flag = false;
    ControlPoint p1, p2;
    p1 = sp.p1;
    p2 = sp.p2;
    
    float [] tt = new float[2];
    float ss;
    PVector p1Old = p1.positionOld.copy();
    PVector p1New = p1.position.copy();
    PVector p2Old = p2.positionOld.copy();
    PVector p2New = p2.position.copy();
    PVector mineOld = cp.positionOld.copy();
    PVector mineNew = cp.position.copy();
    
    p1Old.sub(mineOld);
    p1New.sub(mineNew);
    p2Old.sub(mineOld);
    p2New.sub(mineNew);
    
    PVector a = PVector.sub(new PVector(0,0), p1Old);
    PVector b = PVector.mult(PVector.sub(p1New,p1Old),-1);
    PVector c = PVector.sub(p2Old,p1Old);
    PVector d = PVector.sub(PVector.sub(p2New,p2Old),PVector.sub(p1New,p1Old));
    
    PVector coef2 = b.cross(d);
    PVector coef1 = PVector.add(a.cross(d),b.cross(c));
    PVector coef0 = a.cross(c);
    
    tt = QuadraticRoots( coef2.z, coef1.z, coef0.z );
    
    if (tt[0]>tt[1]) {
      float temp = tt[0];
      tt[0] = tt[1];
      tt[1] = temp;
    }
    
    for (int j=0; j<2 ; j++) {
      if ((tt[j]<0) || (tt[j]>1)) {
        continue;
      }
      else {
        PVector p1Proj = LInterp( p1Old.copy(), p1New.copy(), tt[j]);
        PVector p2Proj = LInterp( p2Old.copy(), p2New.copy(), tt[j]);
        ss = PointProject2Line( p1Proj, p2Proj );
        if ((ss<0) || (ss>1)) {
          continue;
        }
        else {
          PVector cpMove = PVector.sub(cp.position, cp.positionOld);
          PVector p1Move = PVector.sub(p1.position, p1.positionOld);
          PVector p2Move = PVector.sub(p2.position, p2.positionOld);
          float [] tsp = {p1.diameter/(2*p1Move.mag()), p2.diameter/(2*p2Move.mag())};
          this.ResolveCPSpringCollisions( sp, tsp, cp, cp.diameter/(2*cpMove.mag()) );
          Flag = true;
          println("Collision detected");
          println( d );
          delay(2000);
          break;
        }
      }
    }
    return Flag;
  }
  // Resolve control point-spring collisions
  void ResolveCPSpringCollisions( Spring sp, float [] tsp, ControlPoint cp, float tcp ) {
    cp.rewindPosition( tcp );
    sp.p1.rewindPosition( tsp[0] );
    sp.p2.rewindPosition( tsp[1] );
    float otherVelx = 0.5*(sp.p1.velocity.x + sp.p2.velocity.x);
    float otherVely = 0.5*(sp.p1.velocity.y + sp.p2.velocity.y);
    float otherMass = sp.p1.mass + sp.p2.mass;
    
    float Vxi = (cp.velocity.x*(cp.mass-otherMass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*otherVelx;
    float Vyi = (cp.velocity.y*(cp.mass-otherMass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*otherVely;
    float Vxj = (otherVelx*(otherMass-cp.mass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*cp.velocity.x;
    float Vyj = (otherVely*(otherMass-cp.mass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*cp.velocity.y;
    cp.UpdateVelocity( Vxi, Vyi );
    sp.p1.UpdateVelocity( Vxj, Vyj );
    sp.p2.UpdateVelocity( Vxj, Vyj );
  }
  
  // Linear interpolation between two points
  PVector LInterp( PVector start, PVector end, float dt ) {
    PVector out;
    PVector delta = PVector.sub(end, start);
    delta.mult(dt);
    out = PVector.add(start,delta);
    
    return out;
  }
  
  // Project a point onto a line
  float PointProject2Line( PVector start, PVector end ) {
    float s;
    PVector b = PVector.sub(new PVector(0,0), start);
    PVector d = PVector.sub(end, start);
    
    float numer = PVector.dot(b,d);
    float denom = PVector.dot(d,d);
    
    s = numer/denom;
    return s;
  }
  
  // Solve quadratic equation
  float [] QuadraticRoots( float a2, float a1, float a0 ) {
    float [] t = new float[2];
    float dd = sq(a1) - 4*a2*a0;
    if (a2==0) {
      if (a1==0) {
        t[0] = -999;
        t[1] = -999;
      }
      else {
        t[0] = -a0/a1;
        t[1] = t[0];
      }
    }
    else {
      if (dd<0) {
        t[0] = -a1/(2*a2);
        t[1] = t[0];
      }
      else {
        t[0] = (-a1-sqrt(dd))/(2*a2);
        t[1] = (-a1+sqrt(dd))/(2*a2);
      }
    }
    return t;
  }
 
  
  
  // Main Handling method
  void HandleCollisions() {
    // There should be a loop that restarts everytime a collision happened.
    // In this way a collision-free state is ensured by the end of this seemingly 
    // endless loop....
    this.DetectBoundCollision();
    this.DetectCPointCPointCollision();
    this.DetectCPointSpringCollision();
   //this.DetectFastCPCPCollision();
  }
  // Sequential handling of collisions
  void HandleCollisionsSeq() {
    boolean colFlag = true;
    boolean boundFlag = false;
    boolean cpcpFlag = false;
    int iter = 0;
    
    while ((colFlag) || (iter<1000)) {
      iter++;
      for (int i=0; i<Ncp; i++) {
        ControlPoint cp = LocalCPoints.get(i);
        boundFlag = this.DetectBoundCollision( cp );
        if (boundFlag) break;
      }
      for (int i=0; i<Ncp-1; i++) {
        for (int j=i+1; i<Ncp; i++) {
          ControlPoint cpi = LocalCPoints.get(i);
          ControlPoint cpj = LocalCPoints.get(j);
          cpcpFlag = this.DetectCPointCPointCollision( cpi, cpj );
          if (cpcpFlag) break;
        }
      }
      if ((cpcpFlag) || (boundFlag)) colFlag = true; 
    }
  }
  
} // end of class