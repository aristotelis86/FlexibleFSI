class BoundBox {
  Spring spA;
  LineSegment [] lines = new LineSegment[4];
  PVector [] vertices = new PVector[5];
  Window view;
  
  BoundBox( Spring sp1_ ) {
    spA = sp1_;
    this.createAlignedBox();
    view = sp1_.myWindow;
  }
  
  void createAlignedBox() {
    OrthoNormal orth;
    LineSegment line1, line2;
    float xc = 0.5*(spA.p1.position.x + spA.p2.position.x);
    float yc = 0.5*(spA.p1.position.y + spA.p2.position.y);
    
    orth = new OrthoNormal( spA.p1.position, spA.p2.position );
    PVector Normal = new PVector(orth.nx, orth.ny);
    Normal.setMag(spA.p1.diameter/2.);
    
    line1 = new LineSegment( PVector.add(spA.p1.position, Normal), PVector.add(spA.p2.position, Normal));
    line2 = new LineSegment( PVector.sub(spA.p1.position, Normal), PVector.sub(spA.p2.position, Normal));
    
    //mybox = new Body( xc, yc, spA.myWindow );
    //mybox.add(line1.Start.x, line1.Start.y);
    //mybox.add(line1.End.x, line1.End.y);
    //mybox.add(line2.End.x, line2.End.y);
    //mybox.add(line2.Start.x, line2.Start.y);
    //mybox.end();
    lines[0] = line1;
    lines[1] = new LineSegment( new PVector(line1.Start.x, line1.Start.y), new PVector(line2.Start.x, line2.Start.y));
    lines[2] = line2;
    lines[3] = new LineSegment( new PVector(line2.End.x, line2.End.y), new PVector(line1.End.x, line1.End.y));
    vertices[0] = new PVector(xc, yc);
    vertices[1] = new PVector(line1.Start.x, line1.Start.y);
    vertices[2] = new PVector(line1.End.x, line1.End.y);
    vertices[3] = new PVector(line2.Start.x, line2.Start.y);
    vertices[4] = new PVector(line2.End.x, line2.End.y);
  }
  
  void display() {
    stroke(255, 0, 0);
    for (LineSegment ls : lines) {
      line(view.px(ls.Start.x), view.py(ls.Start.y), view.px(ls.End.x), view.py(ls.End.y));
    }
  }
}