class LineSegment {
  PVector Start;
  PVector End;
  OrthoNormal orth;
  
  LineSegment( PVector s, PVector e ) {
    Start = s.copy();
    End = e.copy();
    orth = new OrthoNormal( Start, End );
  }
  
}