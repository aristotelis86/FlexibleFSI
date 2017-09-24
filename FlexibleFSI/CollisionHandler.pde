class CollisionHandler {
  //================= Attributes ====================//
  int Ncp; // number of control points to consider
  ArrayList<ControlPoint> LocalCPoints = new ArrayList<ControlPoint>();
  
  ArrayList<String> BoundCollisionInfo = new ArrayList<String>();
  ArrayList<ControlPoint> BoundCollision = new ArrayList<ControlPoint>();
  
  ArrayList<ControlPoint> CPi = new ArrayList<ControlPoint>();
  ArrayList<ControlPoint> CPj = new ArrayList<ControlPoint>();
  FloatList RewTi = new FloatList();
  FloatList RewTj = new FloatList();
  
  ArrayList<ControlPoint> FastCPi = new ArrayList<ControlPoint>();
  ArrayList<ControlPoint> FastCPj = new ArrayList<ControlPoint>();
  FloatList FastT = new FloatList();
  
  //================= Constructor ====================//
  CollisionHandler( ControlPoint [] cpoints ) {
    Ncp = cpoints.length;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
  }
  CollisionHandler( ControlPoint cpoint ) {
    Ncp = 1;
    LocalCPoints.add(cpoint);
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
            // calculate how much penetration we have....
            // works like that BUT should make a sketch...
            // does not make sense to fully rewind if dij=clearRad
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
 
 // Main Handling method
 void HandleCollisions() {
   this.DetectBoundCollision();
   //this.DetectFastCPCPCollision();
   this.DetectCPointCPointCollision();
   
   this.ResolveNSCollisions();
   //this.ResolveFastCPCPCollisions();
   this.ResolveCPCPCollisions();
 }
  
  
}