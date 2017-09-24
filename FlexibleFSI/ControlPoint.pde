//======================= ControlPoint Class ====================== //

// ------------------- Example Code for Testing ------------------- //

// ---------------------------------------------------------------- //

class ControlPoint {
  //================= Attributes ====================//
  
  PVector position; // current position
  PVector positionOld;
  PVector velocity; // current velocity
  PVector velocityOld;
  PVector acceleration; // current acceleration
  PVector accelerationOld; // current acceleration
  PVector force; // force acting on the point-mass
  PVector impForce; // impact force tracking
  float mass; // mass of the point
  boolean xfixed; // fix the particle at its y-axis
  boolean yfixed; // fix the particle at its x-ayis
  
  Window myWindow; // viewing window
  color c; // for displaying
  float diameter; // for collisions model
  
  //================= Constructor ====================//
  ControlPoint(PVector position_, float m,  float thk, Window myWindow_) {
    position = position_;
    velocity = new PVector(0, 0);
    acceleration = new PVector(0, 0);
    force = new PVector(0, 0);
    mass = m;
    
    xfixed = false;
    yfixed = false;
    
    myWindow = myWindow_;
    c = color(random(1,255), random(1,255), random(1,255));
    diameter = thk;
    
    positionOld = position.copy();
    velocityOld = velocity.copy();
    accelerationOld = acceleration.copy();
  }
  
  ControlPoint(PVector position_, float m, Window myWindow_) {this(position_, m,  1, myWindow_);}
  
  
  //================= Methods =====================//
  
  // Display
  void display() {
    noStroke();
    fill(c);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), diameter*myWindow.x.r, diameter*myWindow.y.r);
  }
  
  // Display impact 
  void impDisplay() {
    noStroke();
    fill(255, 0, 0);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), 0.8*diameter*myWindow.x.r, 0.8*diameter*myWindow.y.r);
  }
  
  // Update the position
  void UpdatePosition(float x, float y) {
    position.x = x;
    position.y = y;
  }
  
  // Update the velocity
  void UpdateVelocity(float x, float y) {
    velocity.x = x;
    velocity.y = y;
  }
  
  // Store old state
  void StoreOld() {
    positionOld = position.copy();
    velocityOld = velocity.copy();
    accelerationOld = acceleration.copy();
  }
  
  // For testing...
  void randomWalk() {
    this.StoreOld();
    float ex, ey, vx, vy, x, y;
    ex = random(1);
    ey = random(1);
    
    if (ex<.05) vx = -1*velocity.x;
    else vx = velocity.x;
    if (ey<.05) vy = -1*velocity.y;
    else vy = velocity.y;
    
    this.UpdateVelocity(vx,vy);
    x = positionOld.x + velocity.x;
    y = positionOld.y + velocity.y;
    this.UpdatePosition(x,y);
  }
  
  // Clear any forces acting on the particle
  void clearForce() { force.mult(0); }
  
  // Accumulate all the forces acting on the particle
  void ApplyForce(PVector FF) { force.add( FF ); }
  
  // Find the acceleration due to forces
  void calculateAcceleration() {
    PVector accel = force.copy();
    acceleration = accel.div(mass);
  }
  
  // Make the particle free of constraints
  void makeFree() {
    xfixed = false;
    yfixed = false;
  }
  void makeFreex() { xfixed = false; }
  void makeFreey() { yfixed = false; }
  
  // Constrain the particle at its location
  void makeFixed() {
    xfixed = true;
    yfixed = true;
  }
  void makeFixedx() { xfixed = true; }
  void makeFixedy() { yfixed = true; }
  
  boolean isFixedx() { return xfixed; }
  boolean isFixedy() { return yfixed; }
  
  // Update methods based on Predictor-Corrector scheme 
  void update( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      StoreOld();
      float x, y, vx, vy;
      
      if (!this.isFixedx()) {
        x = position.x + t*velocity.x;
        vx = velocity.x + t*acceleration.x;
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = position.y + t*velocity.y;
        vy = velocity.y + t*acceleration.y;
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      UpdatePosition( x, y );
      UpdateVelocity( vx, vy );
    }
  }
  
  void update2( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      float x, y, vx, vy;
      if (!this.isFixedx()) {
        x = positionOld.x + .5*t*(velocityOld.x + velocity.x);
        vx = velocityOld.x + .5*t*(accelerationOld.x + acceleration.x);
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = positionOld.y + .5*t*(velocityOld.y + velocity.y);
        vy = velocityOld.y + .5*t*(accelerationOld.y + acceleration.y);
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      UpdatePosition( x, y );
      UpdateVelocity( vx, vy );
    }
  }
  
  // Alternative update methods based on Predictor-Corrector scheme 
  void updateAlt( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      StoreOld();
      float x, y, vx, vy;
      
      if (!this.isFixedx()) {
        x = position.x + t*velocity.x + 0.5*acceleration.x*sq(t);
        vx = velocity.x + t*acceleration.x;
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = position.y + t*velocity.y + 0.5*acceleration.y*sq(t);
        vy = velocity.y + t*acceleration.y;
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      UpdatePosition( x, y );
      UpdateVelocity( vx, vy );
    }
  }
  
  void updateAlt2( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      float vx, vy;
      if (!this.isFixedx()) vx = velocityOld.x + .5*t*(accelerationOld.x + acceleration.x);
      else vx = velocity.x;
      if (!this.isFixedy()) vy = velocityOld.y + .5*t*(accelerationOld.y + acceleration.y);
      else vy = velocity.y;
      UpdateVelocity( vx, vy );
    }
  }
  
  // Export info on the motion of the control points
  void dampInfo( PrintWriter outFile ) {
    outFile.println(position.x+","+position.y+","+velocity.x+","+velocity.y+","+force.x+","+force.y);
  }
  
  // Get the distance between control points
  float distance(ControlPoint other) {
    float d = this.position.dist(other.position);
    return d;
  }
  
  // Update the position given a collision time
  void impUpdate( float dt ) {
    float x, y;
    x = positionOld.x + dt*(position.x - positionOld.x);
    y = positionOld.y + dt*(position.y - positionOld.y);
    UpdatePosition(x,y);
  }
  
  
  
  
  
  
  
  
  
  
  // Resolve CPoint-CPoint collisions event
  void ResolveCPointCPoint ( ControlPoint other, float tt ) {
    float xnewmine, ynewmine, xnewother, ynewother;
    
    xnewmine = positionOld.x + tt*(position.x - positionOld.x);
    ynewmine = positionOld.y + tt*(position.y - positionOld.y);
    xnewother = other.positionOld.x + tt*(other.position.x - other.positionOld.x);
    ynewother = other.positionOld.y + tt*(other.position.y - other.positionOld.y);
    
    float Vxi = (this.velocity.x*(this.mass-other.mass)/(this.mass+other.mass)) + (2*other.mass/(this.mass+other.mass))*other.velocity.x;
    float Vyi = (this.velocity.y*(this.mass-other.mass)/(this.mass+other.mass)) + (2*other.mass/(this.mass+other.mass))*other.velocity.y;
    float Vxj = (other.velocity.x*(other.mass-this.mass)/(this.mass+other.mass)) + (2*other.mass/(this.mass+other.mass))*this.velocity.x;
    float Vyj = (other.velocity.y*(other.mass-this.mass)/(this.mass+other.mass)) + (2*other.mass/(this.mass+other.mass))*this.velocity.y;
    //xnewmine = xnewmine + (1-.6*tt)*Vxi;
    //ynewmine = ynewmine + (1-.6*tt)*Vyi;
    //xnewother = xnewother + (1-.6*tt)*Vxj;
    //ynewother = ynewother + (1-.6*tt)*Vyj;
    this.UpdatePosition( xnewmine, ynewmine );
    other.UpdatePosition( xnewother, ynewother );
    this.UpdateVelocity( Vxi, Vyi );
    other.UpdateVelocity( Vxj, Vyj );
  }
  
  
  
  
  
  
  
  
  
  
  
  //// Apply impact forces to main force variable
  //void applyImpForce() { force.add( impForce ); }
  
  
  
  
  
  
  
  
  //void FastCPointCPointCollision( ControlPoint other ) {
  //  float tol = 1e-6;
  //  float tcol, denom, Rt;
    
  //  float vxmine, vymine, vxother, vyother;
  //  float xoldmine, yoldmine, xoldother, yoldother;
    
  //  float clearRad = (thick/2) + (other.thick/2);
    
  //  vxmine = position.x-positionOld.x;
  //  vymine = position.y-positionOld.y;
  //  vxother = other.position.x-other.positionOld.x;
  //  vyother = other.position.y-other.positionOld.y;
    
  //  xoldmine = positionOld.x; yoldmine = positionOld.y;
  //  xoldother = other.positionOld.x; yoldother = other.positionOld.y;
    
  //  denom = sq(vxmine-vxother)+sq(vymine-vyother);
  //  if (denom<tol) {
  //    CPointCPointCollision( other );
  //  }
  //  else {
  //    tcol = -1*((vxmine-vxother)*(xoldmine-xoldother)+(vymine-vyother)*(yoldmine-yoldother))/denom;
      
  //    Rt = sq((xoldmine+tcol*vxmine) - (xoldother + tcol*vxother)) + sq((yoldmine+tcol*vymine) - (yoldother + tcol*vyother));
      
  //    if (Rt<sq(clearRad)) {
  //      if ((tcol>=0) && (tcol<=1)) {
  //        println("fast cpoint-cpoint collision");
  //        ResolveCPointCPoint( other, tcol );          
  //      }
  //    }
  //  }
  //}
  
  //void CPointSpringCollision( Spring sp ) {
  //  float tol = 1e-7;
  //  float CSDist = -1;
  //  ControlPoint p1, p2;
  //  p1 = sp.p1;
  //  p2 = sp.p2;
    
  //  float clearRad = p1.thick/2 + this.thick/2; 
    
  //  PVector n = PVector.sub(p2.position, p1.position);
  //  PVector cp1 = PVector.sub(p1.position, this.position);
    
  //  float dotProdncp1 = PVector.dot(n,cp1);
    
  //  // Closest point is a
  //  if ( dotProdncp1 > tol ) CSDist = PVector.dot(cp1,cp1);
  //  else {
      
  //    PVector cp2 = PVector.sub(this.position, p2.position);
  //    float dotProdncp2 = PVector.dot(n,cp2);
    
  //    // Closest point is b
  //    if ( dotProdncp2 > tol ) CSDist = PVector.dot(cp2,cp2);
  //    else {
  //      // Closest point is between a and b
  //      PVector epsilon;
  //      float denom = PVector.dot(n,n);
  //      n.mult(dotProdncp1/denom);
  //      epsilon = PVector.sub(cp1,n);
  //      CSDist = PVector.dot(epsilon,epsilon);
  //    }
  //  }
    
  //  if (CSDist >= 0) {
  //    if (sqrt(CSDist)<=clearRad) {
  //      println("clearance="+clearRad);
  //      println("currDist="+sqrt(CSDist));
  //      println("p1="+p1.position.x+" "+p1.position.y);
  //      println("p2="+p2.position.x+" "+p2.position.y);
  //      println("this="+this.position.x+" "+this.position.y);
  //      p1.impDisplay();
  //      p2.impDisplay();
  //      this.impDisplay();
  //      noLoop();
  //    }
  //  }
  //}
  
  
  //void LineSweepsPoint ( Spring sp ) {
  //  ControlPoint p1, p2;
  //  p1 = sp.p1;
  //  p2 = sp.p2;
  //  float [] tt = new float[2];
  //  float ss;
  //  PVector p1Old = p1.positionOld.copy();
  //  PVector p1New = p1.position.copy();
  //  PVector p2Old = p2.positionOld.copy();
  //  PVector p2New = p2.position.copy();
  //  PVector mineOld = this.positionOld.copy();
  //  PVector mineNew = this.position.copy();
    
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
  //        //println("p1="+p1.position.x+" "+p1.position.y);
  //        //println("p2="+p2.position.x+" "+p2.position.y);
  //        //println("this="+this.position.x+" "+this.position.y);
  //        //p1.impDisplay();
  //        //p2.impDisplay();
  //        //this.impDisplay();
  //        //println(tt[j]);
  //        //println(ss);
  //        ResolveCPointSpring( sp, tt[j] );
  //        //noLoop();
  //        break;
  //      }
  //    }
  //  }
  //}
  
  //PVector LInterp( PVector start, PVector end, float dt ) {
  //  PVector out;
    
  //  PVector delta = PVector.sub(end, start);
  //  delta.mult(dt);
  //  out = PVector.add(start,delta);
    
  //  return out;
  //}
  
  //float PointProject2Line( PVector start, PVector end ) {
  //  float s;
  //  PVector b = PVector.sub(new PVector(0,0), start);
  //  PVector d = PVector.sub(end, start);
    
  //  float numer = PVector.dot(b,d);
  //  float denom = PVector.dot(d,d);
    
  //  s = numer/denom;
  //  return s;
  //}
  
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
  
  //void ResolveCPointSpring ( Spring spring, float tt ) {
  //  float xnewmine, ynewmine, xnewother1, ynewother1, xnewother2, ynewother2;
  //  ControlPoint other1 = spring.p1;
  //  ControlPoint other2 = spring.p2;
    
  //  xnewmine = positionOld.x + .5*tt*(position.x - positionOld.x);
  //  ynewmine = positionOld.y + .5*tt*(position.y - positionOld.y);
    
  //  xnewother1 = other1.positionOld.x + .5*tt*(other1.position.x - other1.positionOld.x);
  //  ynewother1 = other1.positionOld.y + .5*tt*(other1.position.y - other1.positionOld.y);
  //  xnewother2 = other2.positionOld.x + .5*tt*(other2.position.x - other2.positionOld.x);
  //  ynewother2 = other2.positionOld.y + .5*tt*(other2.position.y - other2.positionOld.y);
    
  //  float otherVelx = 0.5*(other1.velocity.x + other2.velocity.x);
  //  float otherVely = 0.5*(other1.velocity.y + other2.velocity.y);
  //  float otherMass = other1.mass + other2.mass;
    
    
  //  float Vxi = (this.velocity.x*(this.mass-otherMass)/(this.mass+otherMass)) + (2*otherMass/(this.mass+otherMass))*otherVelx;
  //  float Vyi = (this.velocity.y*(this.mass-otherMass)/(this.mass+otherMass)) + (2*otherMass/(this.mass+otherMass))*otherVely;
  //  float Vxj = (otherVelx*(otherMass-this.mass)/(this.mass+otherMass)) + (2*otherMass/(this.mass+otherMass))*this.velocity.x;
  //  float Vyj = (otherVely*(otherMass-this.mass)/(this.mass+otherMass)) + (2*otherMass/(this.mass+otherMass))*this.velocity.y;
  //  //xnewmine = xnewmine + (1-.6*tt)*Vxi;
  //  //ynewmine = ynewmine + (1-.6*tt)*Vyi;
  //  //xnewother = xnewother + (1-.6*tt)*Vxj;
  //  //ynewother = ynewother + (1-.6*tt)*Vyj;
  //  this.UpdatePosition( xnewmine, ynewmine );
  //  other1.UpdatePosition( xnewother1, ynewother1 );
  //  other2.UpdatePosition( xnewother2, ynewother2 );
  //  this.UpdateVelocity( Vxi, Vyi );
  //  other1.UpdateVelocity( Vxj, Vyj );
  //  other2.UpdateVelocity( Vxj, Vyj );
  //}

} // end of ControlPoint class