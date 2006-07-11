#ifndef FP_H
#define FP_H


struct OptParam
{
  int minNodes_;
  int nodeInterval_;
  double maxTime_;
  OptParam(int minNodes = 100000, int nodeInterval = 1000,double maxTime = 10800):
      minNodes_(minNodes), nodeInterval_(nodeInterval), maxTime_(maxTime)
  {}
}
;

extern OptParam params;


#endif
