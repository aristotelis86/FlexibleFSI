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
  
  ArrayList<String> BoundCollisionInfo = new ArrayList<String>();
  ArrayList<ControlPoint> BoundCollision = new ArrayList<ControlPoint>();
  
  ArrayList<ControlPoint> CPi = new ArrayList<ControlPoint>();
  ArrayList<ControlPoint> CPj = new ArrayList<ControlPoint>();
  FloatList RewTi = new FloatList();
  FloatList RewTj = new FloatList();
  
  ArrayList<ControlPoint> FastCPi = new ArrayList<ControlPoint>();
  ArrayList<ControlPoint> FastCPj = new ArrayList<ControlPoint>();
  FloatList FastT = new FloatList();
  
  ArrayList<ControlPoint> CP = new ArrayList<ControlPoint>();
  ArrayList<Spring> SP = new ArrayList<Spring>();
  FloatList RewTcp = new FloatList();
  FloatList RewTp1 = new FloatList();
  FloatList RewTp2 = new FloatList();
  
  
  //================= Constructor ====================//
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
        BoundCollision.add(cp);
        BoundCollisionInfo.add("West");
      }
      else if (cp.position.x > cp.myWindow.x.inE - cp.diameter/2) {
        BoundCollision.add(cp);
        BoundCollisionInfo.add("East");
      }
      if (cp.position.y < cp.diameter/2) {
        BoundCollision.add(cp);
        BoundCollisionInfo.add("North");
      }
      else if (cp.position.y > cp.myWindow.y.inE - cp.diameter/2) {
        BoundCollision.add(cp);
        BoundCollisionInfo.add("South");
      }
    }
  }
  
  //Resolve North-South Boundary Collisions
  void ResolveNSCollisions() {
    int N = BoundCollision.size();
    if (N>0) {
      for (int i=0; i<N; i++) {
        String str = BoundCollisionInfo.get(i);
        ControlPoint cp = BoundCollision.get(i);
        if ((str.equals("North")==true) || (str.equals("South")==true)) {
          float vx = cp.velocity.x;
          float vy = -cp.velocity.y;
          float x = cp.position.x;
          float y = cp.position.y;
          if (str.equals("North")==true) y = cp.diameter/2;
          else if (str.equals("South")==true) y = cp.myWindow.y.inE-cp.diameter/2;
          cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
        }
        if ((str.equals("West")==true) || (str.equals("East")==true)) {
          float vx = -cp.velocity.x;
          float vy = cp.velocity.y;
          float x = cp.position.x;
          float y = cp.position.y;
          if (str.equals("West")==true) x = cp.diameter/2;
          else if (str.equals("East")==true) x = cp.myWindow.x.inE-cp.diameter/2;
          cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
        }
      }
    }
    BoundCollision = new ArrayList<ControlPoint>();
    BoundCollisionInfo = new ArrayList<String>();
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
            CPi.add(cpi);
            CPj.add(cpj);
            float penet = 0.5*(cpi.diameter + cpj.diameter) - dij;
            float rewindi = penet/(0.5*cpi.diameter);
            float rewindj = penet/(0.5*cpj.diameter);
            RewTi.append(rewindi);
            RewTj.append(rewindj);
            // Must include a detection method for fast moving cpoints
            // Cheking if their relative velocity is greater than clearRad 
          }
        }
      }
    }
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
  
  // Resolve CPoint-CPoint collisions
  void ResolveCPCPCollisions() {
    int N = CPi.size();
    for (int ij=0; ij<N; ij++) {
      ControlPoint cpi = CPi.get(ij);
      ControlPoint cpj = CPj.get(ij);
      float rewTi = RewTi.get(ij);
      float rewTj = RewTj.get(ij);
      cpi.rewindPosition(rewTi);
      cpj.rewindPosition(rewTj);
      float Vxi = (cpi.velocity.x*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.x;
      float Vyi = (cpi.velocity.y*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.y;
      float Vxj = (cpj.velocity.x*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.x;
      float Vyj = (cpj.velocity.y*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.y;
      cpi.UpdateVelocity( Vxi, Vyi );
      cpj.UpdateVelocity( Vxj, Vyj );
    }
    CPi = new ArrayList<ControlPoint>();
    CPj = new ArrayList<ControlPoint>();
    RewTi = new FloatList();
    RewTj = new FloatList();
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
  
  void LineSweepsPoint( Spring sp, ControlPoint cp ) {
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
          CP.add(cp);
          SP.add(sp);
          RewTcp.append(cp.diameter/(2*cpMove.mag()));
          RewTp1.append(p1.diameter/(2*p1Move.mag()));
          RewTp2.append(p2.diameter/(2*p2Move.mag()));
          println(cp.diameter/(2*cpMove.mag()));
          println(p1.diameter/(2*p1Move.mag()));
          println(p2.diameter/(2*p2Move.mag()));
          break;
        }
      }
    }
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
  
  // Resolve control point-spring collisions
  // Problem resolving simultaneous collisions must not correct multiple times....
  // Maybe check for parallel spring...?
  void ResolveCPSpringCollisions() {
    int N = CP.size();
    if (N>0) {
      for (int i=0; i<N; i++) {
        ControlPoint cp = CP.get(i);
        Spring sp = SP.get(i);
        println("cp before "+cp.position.x+", "+cp.position.y);
        println("p1 before "+sp.p1.position.x+", "+sp.p1.position.y);
        println("p2 before "+sp.p2.position.x+", "+sp.p2.position.y);
        cp.rewindPosition( RewTcp.get(i) );
        //sp.p1.rewindPosition( RewTp1.get(i) );
        sp.p2.rewindPosition( RewTp2.get(i) );
        println("cp after "+cp.position.x+", "+cp.position.y);
        println("p1 after "+sp.p1.position.x+", "+sp.p1.position.y);
        println("p2 after "+sp.p2.position.x+", "+sp.p2.position.y);
        noLoop();
      }
      CP = new ArrayList<ControlPoint>();
      SP = new ArrayList<Spring>();
      RewTcp = new FloatList();
      RewTp1 = new FloatList();
      RewTp2 = new FloatList();
    }
  }
  
 
 // Main Handling method
 void HandleCollisions() {
   this.DetectBoundCollision();
   //this.DetectFastCPCPCollision();
   this.DetectCPointCPointCollision();
   this.DetectCPointSpringCollision();
   
   this.ResolveNSCollisions();
   //this.ResolveFastCPCPCollisions();
   this.ResolveCPCPCollisions();
   this.ResolveCPSpringCollisions();
 }
  
  
}