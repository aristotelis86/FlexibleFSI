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
  CollisionHandler( FlexibleSheet [] sheet ) {
    int Nfs = sheet.length;
    
    Ncp = 0;
    Nsp = 0;
    for (int ij=0; ij<Nfs; ij++) {
      for (int i=0; i<sheet[ij].numOfpoints; i++) {
        LocalCPoints.add(sheet[ij].cpoints[i]);
        Ncp++;
      }
      for (int i=0; i<sheet[ij].numOfsprings; i++) {
        LocalSprings.add(sheet[ij].springs[i]);
        Nsp++;
      }
    }
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
  void DetectBoundCollision( ControlPoint cp ) {
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
  // Resolve boundary collisions
  void ResolveBoundCollisions( String bound, ControlPoint cp ){
    float e = 0.4;
    if ((bound.equals("North")==true) || (bound.equals("South")==true)) {
      float vx = cp.velocity.x;
      float vy = -e*cp.velocity.y;
      float x = cp.position.x;
      float y = cp.position.y;
      if (bound.equals("North")==true) y = cp.diameter/2;
      else if (bound.equals("South")==true) y = cp.myWindow.y.inE-cp.diameter/2;
      cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
    }
    if ((bound.equals("West")==true) || (bound.equals("East")==true)) {
      float vx = -e*cp.velocity.x;
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
  
  // Detect spring-spring intersection
  void DetectSpringSpringCollision() {
    if (Nsp>1) {
      for (int i=0; i<Nsp-1; i++) {
        for (int j=i+1; j<Nsp; j++) {
          Spring spi = LocalSprings.get(i);
          Spring spj = LocalSprings.get(j);
          ControlPoint cpi1 = spi.p1;
          ControlPoint cpj1 = spj.p1;
          ControlPoint cpi2 = spi.p2;
          ControlPoint cpj2 = spj.p2;
          if ((((cpi1!=cpj1) && (cpi1!=cpj2))) &&  (((cpi2!=cpj1) && (cpi2!=cpj2)))) {
            BoundBox Boxi, Boxj;
            Boxi = new BoundBox( spi );
            Boxj = new BoundBox( spj );
            //Boxi.displayOBB();
            //Boxj.displayOBB();
            boolean flag = this.doBoundBoxesOverlap( Boxi, Boxj );
            if (flag) {
              this.ResolveSpringSpringCollisions( spi, spj );
              //noLoop();
            }
          }
        }
      }
    }
  }
  void ResolveSpringSpringCollisions( Spring sp1, Spring sp2 ) {
    
      float sp1Velx = 0.5*(sp1.p1.velocity.x + sp1.p2.velocity.x);
      float sp1Vely = 0.5*(sp1.p1.velocity.y + sp1.p2.velocity.y);
      float sp1Mass = sp1.p1.mass + sp1.p2.mass;
      
      float sp2Velx = 0.5*(sp2.p1.velocity.x + sp2.p2.velocity.x);
      float sp2Vely = 0.5*(sp2.p1.velocity.y + sp2.p2.velocity.y);
      float sp2Mass = sp2.p1.mass + sp2.p2.mass;
      
      float Vx1 = (sp1Velx*(sp1Mass-sp2Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp2Velx;
      float Vy1 = (sp1Vely*(sp1Mass-sp2Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp2Vely;
      float Vx2 = (sp2Velx*(sp2Mass-sp1Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp1Velx;
      float Vy2 = (sp2Vely*(sp2Mass-sp1Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp1Vely;
      sp1.p1.UpdateVelocity( Vx1, Vy1 );
      sp1.p2.UpdateVelocity( Vx1, Vy1 );
      sp2.p1.UpdateVelocity( Vx2, Vy2 );
      sp2.p2.UpdateVelocity( Vx2, Vy2 );
  }
  
  boolean doBoundBoxesOverlap( BoundBox b1, BoundBox b2 ) {
    Boolean overFlag = false;
    PVector [] CheckAxes = new PVector[4];
    Boolean [] AxesColFlag = new Boolean[4];
    CheckAxes[0] = new PVector(b1.lines[0].orth.nx, b1.lines[0].orth.ny);
    CheckAxes[1] = new PVector(b1.lines[1].orth.nx, b1.lines[1].orth.ny);
    CheckAxes[2] = new PVector(b2.lines[0].orth.nx, b2.lines[0].orth.ny);
    CheckAxes[3] = new PVector(b2.lines[1].orth.nx, b2.lines[1].orth.ny);
    for (Boolean ax : AxesColFlag) {
      ax = false;
    }
    
    axesLoop : for (int i=0; i<4; i++) {
      float [] ProjB1 = new float[5];
      float [] ProjB2 = new float[5];
      for (int j=0; j<5; j++) {
        ProjB1[j] = PVector.dot(CheckAxes[i],b1.vertices[j]);
        ProjB2[j] = PVector.dot(CheckAxes[i],b2.vertices[j]);
      }
      float minProjB1 = ProjB1[1];
      float maxProjB1 = ProjB1[1];
      float minProjB2 = ProjB2[1];
      float maxProjB2 = ProjB2[1];
      
      for (int j=2; j<5; j++) {
        if (ProjB1[j]<minProjB1) minProjB1 = ProjB1[j];
        if (ProjB1[j]>maxProjB1) maxProjB1 = ProjB1[j];
        if (ProjB2[j]<minProjB2) minProjB2 = ProjB2[j];
        if (ProjB2[j]>maxProjB2) maxProjB2 = ProjB2[j];
      }
      if ((maxProjB2 < minProjB1) || (maxProjB1 < minProjB2)) {
        AxesColFlag[i] = false;
        break axesLoop;
      }
      else AxesColFlag[i] = true;
    }
    if (AxesColFlag[0] && AxesColFlag[1] && AxesColFlag[2] && AxesColFlag[3]) {
      overFlag = true;
    }
    return overFlag;
  }
  
  
  
  //boolean doBoundingBoxesIntersect( BoundBox a, BoundBox b ) {
  //  return a.A[0].y <= b.A[1].y 
  //        && a.A[1].y >= b.A[0].y 
  //        && a.A[1].x >= b.A[0].x 
  //        && a.A[0].x <= b.A[1].x;
  //}
  //boolean isPointOnLine( LineSegment a, PVector b ) {
  //  float tol = 1e-7;
  //  LineSegment aTmp = new LineSegment(new PVector(0, 0), new PVector(
  //          a.End.x - a.Start.x, a.End.y - a.Start.y));
  //  PVector bTmp = new PVector(b.x - a.Start.x, b.y - a.Start.y);
  //  PVector r = aTmp.End.cross( bTmp );
  //  return abs(r.z) < tol;
  //}
  //boolean isPointRightOfLine( LineSegment a, PVector b ) {
  //  LineSegment aTmp = new LineSegment(new PVector(0, 0), new PVector(
  //          a.End.x - a.Start.x, a.End.y - a.Start.y));
  //  PVector bTmp = new PVector(b.x - a.Start.x, b.y - a.Start.y);
  //  PVector r = aTmp.End.cross( bTmp ); 
  //  return r.z < 0;
  //}
  //boolean lineSegmentTouchesOrCrossesLine( LineSegment a, LineSegment b ) {
  //  return isPointOnLine(a, b.Start)
  //          || isPointOnLine(a, b.End)
  //          || (isPointRightOfLine(a, b.Start) ^ isPointRightOfLine(a, b.End));
  //}
  //boolean doLinesIntersect( Spring sp1, Spring sp2 ) {
  //  BoundBox box1 = new BoundBox( sp1 );
  //  BoundBox box2 = new BoundBox( sp2 );
  //  //box1.display();
  //  //box2.display();
  //  box1.mybox.display(#00FF00);
  //  box2.mybox.display(#00FF00);
  //  LineSegment s1 = new LineSegment( sp1.p1.position, sp1.p2.position );
  //  LineSegment s2 = new LineSegment( sp2.p1.position, sp2.p2.position );
    
  //  return doBoundingBoxesIntersect(box1, box2);
  //          //&& lineSegmentTouchesOrCrossesLine(s1, s2)
  //          //&& lineSegmentTouchesOrCrossesLine(s2, s1);
  //}
  
  
  
  
  //// Maybe check for parallel spring...?
  //// Still some false-positive checks occuring
  //// Detect point-spring collisions
  //void  DetectCPointSpringCollision() {
  //  if ((Nsp>0) && (Ncp>2)) {
  //    for (int icp=0; icp<Ncp; icp++) {
  //      ControlPoint cp = LocalCPoints.get(icp);
  //      for (int isp=0; isp<Nsp; isp++) {
  //        Spring sp = LocalSprings.get(isp);
  //        ControlPoint p1, p2;
  //        p1 = sp.p1;
  //        p2 = sp.p2;
  //        if ((cp!=p1) && (cp!=p2)) {
  //          this.LineSweepsPoint( sp, cp );
  //        }
  //      }
  //    }
  //  }
  //}
  // Detect point-spring collisions
  //boolean  DetectCPointSpringCollision( Spring sp, ControlPoint cp ) {
  //  boolean Flag = false;
  //  ControlPoint p1, p2;
  //  p1 = sp.p1;
  //  p2 = sp.p2;
  //  if ((cp!=p1) && (cp!=p2)) {
  //    Flag = LineSweepsPoint( sp, cp );
  //  }
  //  return Flag;
  //}
  // Check if a line sweeps a point
  //void LineSweepsPoint( Spring sp, ControlPoint cp ) {
  //  ControlPoint p1, p2;
  //  p1 = sp.p1;
  //  p2 = sp.p2;
    
  //  float [] tt = new float[2];
  //  float ss;
  //  PVector p1Old = p1.positionOld.copy();
  //  PVector p1New = p1.position.copy();
  //  PVector p2Old = p2.positionOld.copy();
  //  PVector p2New = p2.position.copy();
  //  PVector mineOld = cp.positionOld.copy();
  //  PVector mineNew = cp.position.copy();
    
  //  p1Old.sub(mineOld);
  //  p1New.sub(mineNew);
  //  p2Old.sub(mineOld);
  //  p2New.sub(mineNew);
    
  //  PVector a = PVector.sub(new PVector(0,0), p1Old);
  //  PVector b = PVector.mult(PVector.sub(p1New,p1Old),-1);
  //  PVector c = PVector.sub(p2Old,p1Old);
  //  PVector d = PVector.sub(PVector.sub(p2New,p2Old),PVector.sub(p1New,p1Old));
    
  //  PVector coef2 = b.cross(d);
  //  PVector coef1 = PVector.add(a.cross(d),b.cross(c));
  //  PVector coef0 = a.cross(c);
    
  //  tt = QuadraticRoots( coef2.z, coef1.z, coef0.z );
    
  //  if (tt[0]>tt[1]) {
  //    float temp = tt[0];
  //    tt[0] = tt[1];
  //    tt[1] = temp;
  //  }
    
  //  for (int j=0; j<2 ; j++) {
  //    if ((tt[j]<0) || (tt[j]>1)) {
  //      continue;
  //    }
  //    else {
  //      PVector p1Proj = LInterp( p1Old.copy(), p1New.copy(), tt[j]);
  //      PVector p2Proj = LInterp( p2Old.copy(), p2New.copy(), tt[j]);
  //      ss = PointProject2Line( p1Proj, p2Proj );
  //      if ((ss<0) || (ss>1)) {
  //        continue;
  //      }
  //      else {
  //        PVector cpMove = PVector.sub(cp.position, cp.positionOld);
  //        PVector p1Move = PVector.sub(p1.position, p1.positionOld);
  //        PVector p2Move = PVector.sub(p2.position, p2.positionOld);
  //        float [] tsp = {p1.diameter/(2*p1Move.mag()), p2.diameter/(2*p2Move.mag())};
  //        //float [] tsp = {tt[j], tt[j]};
  //        this.ResolveCPSpringCollisions( sp, tsp, cp, cp.diameter/(2*cpMove.mag()) );
  //        //this.ResolveCPSpringCollisions( sp, tsp, cp, tt[j] );
  //        println("Collision detected");
  //        println( d );
  //        //noLoop();
  //        break;
  //      }
  //    }
  //  }
  //}
  //// Resolve control point-spring collisions
  //void ResolveCPSpringCollisions( Spring sp, float [] tsp, ControlPoint cp, float tcp ) {
  //  cp.rewindPosition( tcp );
  //  sp.p1.rewindPosition( tsp[0] );
  //  sp.p2.rewindPosition( tsp[1] );
  //  float otherVelx = 0.5*(sp.p1.velocity.x + sp.p2.velocity.x);
  //  float otherVely = 0.5*(sp.p1.velocity.y + sp.p2.velocity.y);
  //  float otherMass = sp.p1.mass + sp.p2.mass;
    
  //  float Vxi = (cp.velocity.x*(cp.mass-otherMass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*otherVelx;
  //  float Vyi = (cp.velocity.y*(cp.mass-otherMass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*otherVely;
  //  float Vxj = (otherVelx*(otherMass-cp.mass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*cp.velocity.x;
  //  float Vyj = (otherVely*(otherMass-cp.mass)/(cp.mass+otherMass)) + (2*otherMass/(cp.mass+otherMass))*cp.velocity.y;
  //  cp.UpdateVelocity( Vxi, Vyi );
  //  sp.p1.UpdateVelocity( Vxj, Vyj );
  //  sp.p2.UpdateVelocity( Vxj, Vyj );
  //}
  
  //// Linear interpolation between two points
  //PVector LInterp( PVector start, PVector end, float dt ) {
  //  PVector out;
  //  PVector delta = PVector.sub(end, start);
  //  delta.mult(dt);
  //  out = PVector.add(start,delta);
    
  //  return out;
  //}
  
  //// Project a point onto a line
  //float PointProject2Line( PVector start, PVector end ) {
  //  float s;
  //  PVector b = PVector.sub(new PVector(0,0), start);
  //  PVector d = PVector.sub(end, start);
    
  //  float numer = PVector.dot(b,d);
  //  float denom = PVector.dot(d,d);
    
  //  s = numer/denom;
  //  return s;
  //}
  
  //// Solve quadratic equation
  //float [] QuadraticRoots( float a2, float a1, float a0 ) {
  //  float [] t = new float[2];
  //  float dd = sq(a1) - 4*a2*a0;
  //  if (a2==0) {
  //    if (a1==0) {
  //      t[0] = -999;
  //      t[1] = -999;
  //    }
  //    else {
  //      t[0] = -a0/a1;
  //      t[1] = t[0];
  //    }
  //  }
  //  else {
  //    if (dd<0) {
  //      t[0] = -a1/(2*a2);
  //      t[1] = t[0];
  //    }
  //    else {
  //      t[0] = (-a1-sqrt(dd))/(2*a2);
  //      t[1] = (-a1+sqrt(dd))/(2*a2);
  //    }
  //  }
  //  return t;
  //}
  
  //// Detect Fast moving CPoint-CPoint collisions
  //void DetectFastCPCPCollision() {
  //  int N = LocalCPoints.size();
  //  if (N>1) {
  //    for (int i=0; i<N-1; i++) {
  //      for (int j=i+1; j<N; j++) {
  //        ControlPoint cpi = LocalCPoints.get(i);
  //        ControlPoint cpj = LocalCPoints.get(j);
  //        float tol = 1e-7;
  //        float tcol, denom, Rt;
  //        float vxi, vyi, vxj, vyj;
  //        float xoldi, yoldi, xoldj, yoldj;
  //        vxi = cpi.position.x-cpi.positionOld.x;
  //        vyi = cpi.position.y-cpi.positionOld.y;
  //        vxj = cpj.position.x-cpj.positionOld.x;
  //        vyj = cpj.position.y-cpj.positionOld.y;
  //        xoldi = cpi.positionOld.x; yoldi = cpi.positionOld.y;
  //        xoldj = cpj.positionOld.x; yoldj = cpj.positionOld.y;
  //        denom = sq(vxi-vxj)+sq(vyi-vyj);
  //        if (denom<tol) {
  //          // Moving in parallel or not moving at all....
  //        }
  //        else {
  //          tcol = -1*((vxi-vxj)*(xoldi-xoldj)+(vyi-vyj)*(yoldi-yoldj))/denom;
  //          Rt = sq((xoldi+tcol*vxi) - (xoldj + tcol*vxj)) + sq((yoldi+tcol*vyi) - (yoldj + tcol*vyj));
  //          float clearRad = (cpi.diameter + cpj.diameter)/2.;
  //          if (Rt<sq(clearRad)) {
  //            if ((tcol>=0) && (tcol<=1)) {
  //              // Fast collisions happened!
  //              FastCPi.add(cpi);
  //              FastCPj.add(cpj);
  //              FastT.append(tcol);
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}
  
  //// Resolve Fast moving CPoint-CPoint collisions
  //void ResolveFastCPCPCollisions() {
  //  int N = FastCPi.size();
  //  for (int ij=0; ij<N; ij++) {
  //    ControlPoint cpi = FastCPi.get(ij);
  //    ControlPoint cpj = FastCPj.get(ij);
  //    float tcol = FastT.get(ij);
      
  //    float Vxi = (cpi.velocity.x*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.x;
  //    float Vyi = (cpi.velocity.y*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.y;
  //    float Vxj = (cpj.velocity.x*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.x;
  //    float Vyj = (cpj.velocity.y*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.y;
  //    cpi.UpdateVelocity( Vxi, Vyi );
  //    cpj.UpdateVelocity( Vxj, Vyj );
  //    //cpi.impUpdate(tcol);
  //    //cpj.impUpdate(tcol);
  //    delay(3000);
  //  }
  //  FastCPi = new ArrayList<ControlPoint>();
  //  FastCPj = new ArrayList<ControlPoint>();
  //  FastT = new FloatList();
  //}
 
  
  
  // Main Handling method
  void HandleCollisions() {
    // There should be a loop that restarts everytime a collision happened.
    // In this way a collision-free state is ensured by the end of this seemingly 
    // endless loop....
    this.DetectBoundCollision();
    this.DetectCPointCPointCollision();
    this.DetectSpringSpringCollision();
  }
  //// Sequential handling of collisions
  //void HandleCollisionsSeq() {
  //  boolean colFlag = true;
  //  boolean boundFlag = false;
  //  boolean cpcpFlag = false;
  //  int iter = 0;
    
  //  for (int i=0; i<Ncp; i++) {
  //    ControlPoint cp = LocalCPoints.get(i);
  //    boundFlag = this.DetectBoundCollision( cp );
  //  }
  //  for (int i=0; i<Ncp-1; i++) {
  //    for (int j=i+1; i<Ncp; i++) {
  //      ControlPoint cpi = LocalCPoints.get(i);
  //      ControlPoint cpj = LocalCPoints.get(j);
  //      cpcpFlag = this.DetectCPointCPointCollision( cpi, cpj );
  //    }
  //  }
  //}
  
} // end of class