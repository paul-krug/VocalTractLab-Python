// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#include "GesturalScore.h"
#include "Dsp.h"
#include "Constants.h"
#include "Sampa.h"

#include <cstdio>
#include <cstdlib>

// Reference frequency for the conversion between Hz and st.
const double GesturalScore::REFERENCE_FREQUENCY = 1.0;
const double GesturalScore::DEFAULT_TIME_CONSTANT_S = 0.012;

static const double TARGET_VELIC_OPENING = 0.5;

// ****************************************************************************
/// Constructor. Init the variables.
// ****************************************************************************

GestureSequence::GestureSequence()
{
  name = "";
  abbr = "";
  unit = "";
  minValue = 0.0;
  maxValue = 1.0;
  minSlope = 0.0;
  maxSlope = 0.0;
  minTau_s = 0.0;
  maxTau_s = 0.0;
  nominalValues = false;

  gesture.clear();
}


// ****************************************************************************
/// Sets the common properties for the gestures in this sequence.
// ****************************************************************************

void GestureSequence::init(const string &name, const string &abbr, const string &unit, 
  double minValue, double maxValue, double minSlope, double maxSlope, 
  double minTau_s, double maxTau_s, bool nominalValues)
{
  this->name = name;
  this->abbr = abbr;
  this->unit = unit;
  this->minValue = minValue;
  this->maxValue = maxValue;
  this->minSlope = minSlope;
  this->maxSlope = maxSlope;
  this->minTau_s = minTau_s;
  this->maxTau_s = maxTau_s;
  this->nominalValues = nominalValues;
}


// ****************************************************************************
/// Clear the sequence.
// ****************************************************************************

void GestureSequence::clear()
{
  gesture.clear();
}


// ****************************************************************************
/// Returns the number of gestures in this sequence.
// ****************************************************************************

int GestureSequence::numGestures()
{
  return (int)gesture.size();
}


// ****************************************************************************
/// Returns a pointer to the gesture with the given index or NULL, if the index
/// is out of range.
// ****************************************************************************

Gesture *GestureSequence::getGesture(int index)
{
  if (isValidIndex(index))
  {
    return &gesture[index];
  }
  else
  {
    return NULL;
  }
}


// ****************************************************************************
/// Returns true, if the index is that of an existing gesture.
// ****************************************************************************

bool GestureSequence::isValidIndex(int index)
{
  return ((index >= 0) && (index < (int)gesture.size()));
}


// ****************************************************************************
/// Returns the start time of the gesture with the given index, or 0, if the
/// index is out of range.
// ****************************************************************************

double GestureSequence::getGestureBegin_s(int index)
{
  double t_s = 0.0;
  
  if (isValidIndex(index))
  {
    int i;
    for (i=0; i < index; i++)
    {
      t_s+= gesture[i].duration_s;      
    }
  }

  return t_s;
}


// ****************************************************************************
/// Returns the end time of the gesture with the given index, or 0, if the
/// index is out of range.
// ****************************************************************************

double GestureSequence::getGestureEnd_s(int index)
{
  double t_s = 0.0;
  
  if (isValidIndex(index))
  {
    int i;
    for (i=0; i <= index; i++)
    {
      t_s+= gesture[i].duration_s;      
    }
  }

  return t_s;
}


// ****************************************************************************
/// Returns the index of the gesture at the position pos_s or -1, if there is
/// no gesture.
// ****************************************************************************

int GestureSequence::getIndexAt(double pos_s)
{
  int i, k = -1;
  double start_s = 0.0;

  for (i=0; (i < (int)gesture.size()) && (k == -1); i++)
  {
    if ((pos_s >= start_s) && (pos_s < start_s + gesture[i].duration_s)) 
    { 
      k = i; 
    }
    start_s+= gesture[i].duration_s;
  }

  return k;
}


// ****************************************************************************
/// Returns the duration of this sequence.
// ****************************************************************************

double GestureSequence::getDuration_s()
{
  int i;
  double duration_s = 0.0;

  for (i=0; i < (int)gesture.size(); i++)
  {
    duration_s+= gesture[i].duration_s;
  }

  return duration_s;
}


// ****************************************************************************
/// Adds the given gesture to the end of the sequence.
// ****************************************************************************

void GestureSequence::appendGesture(Gesture &g)
{
  gesture.push_back(g);
}


// ****************************************************************************
/// Inserts a new gesture g at the position index into the sequence.
// ****************************************************************************

void GestureSequence::insertGesture(Gesture &g, int index)
{
  if (isValidIndex(index) == false) 
  { 
    return; 
  }
  vector<Gesture>::iterator iter = gesture.begin();
  advance(iter, index);
  gesture.insert(iter, g);
}


// ****************************************************************************
/// Deletes the gesture at the position index.
// ****************************************************************************

void GestureSequence::deleteGesture(int index)
{
  if (isValidIndex(index) == false) 
  { 
    return; 
  }
  vector<Gesture>::iterator iter = gesture.begin();
  advance(iter, index);
  gesture.erase(iter);
}


// ****************************************************************************
/// Puts the gesture g at the position startPos. Underlying gestures are
/// overwritten or skipped accordingly.
// ****************************************************************************

void GestureSequence::putGesture(Gesture &g, double startPos_s)
{
  const double MIN_GESTURE_DURATION_S = (double)MIN_GESTURE_DURATION_MS / 1000.0;

  // Range check.
  if (g.duration_s < MIN_GESTURE_DURATION_S)
  {
    g.duration_s = MIN_GESTURE_DURATION_S;
  }

  // End position of the new gesture.
  double endPos_s = startPos_s + g.duration_s;

  int i;
  int leftGestureIndex = -1;
  int rightGestureIndex = -1;
  double leftGestureStart_s = 0.0;
  double rightGestureStart_s = 0.0;
  double pos_s = 0.0;

  for (i = 0; i < (int)gesture.size(); i++)
  {
    if ((startPos_s >= pos_s) && (startPos_s <= pos_s + gesture[i].duration_s))
    {
      leftGestureIndex = i;
      leftGestureStart_s = pos_s;
    }

    if ((endPos_s >= pos_s) && (endPos_s <= pos_s + gesture[i].duration_s))
    {
      rightGestureIndex = i;
      rightGestureStart_s = pos_s;
    }

    pos_s += gesture[i].duration_s;
  }

  double totalDuration_s = pos_s;

  // ****************************************************************
  // The start of the new gesture is after the last existing gesture.
  // ****************************************************************

  if ((leftGestureIndex == -1) && (rightGestureIndex == -1))
  {
    double emptyDuration_s = startPos_s - totalDuration_s;
    if (emptyDuration_s > MIN_GESTURE_DURATION_S)
    {
      Gesture emptyGesture;
      emptyGesture.duration_s = emptyDuration_s;
      emptyGesture.dVal = 0.0;
      emptyGesture.slope = 0.0;
      emptyGesture.sVal = "";
      emptyGesture.neutral = true;
      emptyGesture.tau_s = GesturalScore::DEFAULT_TIME_CONSTANT_S;

      appendGesture(emptyGesture);
    }

    // Now append the intended new gesture.
    appendGesture(g);
  }
  else

  // ****************************************************************
  // The new gesture overlaps at least partly with existing gestures.
  // ****************************************************************

  if (leftGestureIndex != -1)     // Should here always be the case.
  {
    // Create the new vector of gestures.
    vector<Gesture> newGestures;

    // Copy all gestures of the original sequence up to leftGestureIndex-1
    for (i = 0; i < leftGestureIndex; i++)
    {
      newGestures.push_back(gesture[i]);
    }

    // Append the beginning of the original gesture leftGestureIndex

    Gesture gesturePart = gesture[leftGestureIndex];
    gesturePart.duration_s = startPos_s - leftGestureStart_s;

    if (gesturePart.duration_s >= MIN_GESTURE_DURATION_S)
    {
      newGestures.push_back(gesturePart);
    }
    else
    {
      // Add the little duration of the of the deleted gesture part
      // to the duration of the new gesture.
      g.duration_s += gesturePart.duration_s;
    }

    // Append the new gesture.
    newGestures.push_back(g);

    if (rightGestureIndex != -1)
    {
      // Append the remaining part of the rightGestureIndex
      gesturePart = gesture[rightGestureIndex];
      gesturePart.duration_s = rightGestureStart_s + gesture[rightGestureIndex].duration_s - endPos_s;
      if (gesturePart.duration_s < MIN_GESTURE_DURATION_S)
      {
        gesturePart.duration_s = MIN_GESTURE_DURATION_S;
      }
      newGestures.push_back(gesturePart);

      // Append all original gestures from rightGestureIndex+1 on.
      for (i = rightGestureIndex + 1; i < (int)gesture.size(); i++)
      {
        newGestures.push_back(gesture[i]);
      }
    }

    // Substitute the original gesture sequence with the new sequence.
    gesture.clear();
    gesture = newGestures;
  }
}


// ****************************************************************************
/// Writes the data of this gesture sequence into the given ASCII stream in
/// XML format.
// ****************************************************************************

bool GestureSequence::writeToXml(ostream &os, int indent)
{
  char st[1024];
  string value;
  int i;

  // ****************************************************************
  // Output the start-tag of this element.
  // ****************************************************************

  sprintf(st, "<gesture_sequence type=\"%s\" unit=\"%s\">", abbr.c_str(), unit.c_str());
  os << string(indent, ' ') << st << endl;

  // ****************************************************************
  // Output one individual element for each gesture.
  // ****************************************************************

  indent+= 2;

  for (i=0; i < (int)gesture.size(); i++)
  {
    if (nominalValues)
    {
      value = gesture[i].sVal;
    }
    else
    {
      sprintf(st, "%f", gesture[i].dVal);
      value = string(st);
    }

    sprintf(st, "<gesture value=\"%s\" slope=\"%f\" duration_s=\"%f\" "
      "time_constant_s=\"%f\" neutral=\"%d\" />",
      value.c_str(),
      gesture[i].slope,
      gesture[i].duration_s,
      gesture[i].tau_s,
      (int)gesture[i].neutral);

    os << string(indent, ' ') << st << endl;
  }

  indent-= 2;

  // ****************************************************************
  // End-tag of the sequence element.
  // ****************************************************************
  
  os << string(indent, ' ') << "</gesture_sequence>" << endl;
  
  return true;
}


// ****************************************************************************
/// Reads the data for this gesture sequence from the given XML tree.
// ****************************************************************************

bool GestureSequence::readFromXml(XmlNode &node, bool &allValuesInRange)
{
  int i;
  int numChilds = node.numChildElements("gesture");
  XmlNode *gestureNode;
  Gesture g;
  string value;

  // Clear all existing gestures in this sequence.
  clear();
  // Assume the best:
  allValuesInRange = true;

  // ****************************************************************
  // Run through all gestures in the sequence.
  // ****************************************************************

  for (i=0; i < numChilds; i++)
  {
    gestureNode = node.getChildElement("gesture", i);
    
    // Extract the gesture properties.
    
    value = gestureNode->getAttributeString("value");
    if (nominalValues)
    {
      g.sVal = value;
      g.dVal = 0.0;
    }
    else
    {
      g.sVal = "";
      g.dVal = atof(value.c_str());
    
      // Check the value range
      if (g.dVal < minValue)
      {
        g.dVal = minValue;
        printf("Gesture value was too low and has been corrected.\n");
        allValuesInRange = false;
      }
      if (g.dVal > maxValue)
      {
        g.dVal = maxValue;
        printf("Gesture value was too high and has been corrected.\n");
        allValuesInRange = false;
      }
    }

    // **************************************************************

    g.slope      = gestureNode->getAttributeDouble("slope");
    g.duration_s = gestureNode->getAttributeDouble("duration_s");
    g.tau_s      = gestureNode->getAttributeDouble("time_constant_s");

    // Check the ranges.
    
    if (g.slope < minSlope)
    {
      g.slope = minSlope;
      printf("Gesture slope was too low and has been corrected.\n");
      allValuesInRange = false;
    }
    if (g.slope > maxSlope)
    {
      g.slope = maxSlope;
      printf("Gesture slope was too high and has been corrected.\n");
      allValuesInRange = false;
    }

    if (g.duration_s < (double)MIN_GESTURE_DURATION_MS / 1000.0)
    {
      g.duration_s = (double)MIN_GESTURE_DURATION_MS / 1000.0;
      printf("Gesture duration was too low and has been corrected.\n");
      allValuesInRange = false;
    }
    if (g.duration_s > (double)MAX_GESTURE_DURATION_MS / 1000.0)
    {
      g.duration_s = (double)MAX_GESTURE_DURATION_MS / 1000.0;
      printf("Gesture duration was too high and has been corrected.\n");
      allValuesInRange = false;
    }

    if (g.tau_s < minTau_s)
    {
      g.tau_s = minTau_s;
      printf("Gesture time constant was too low and has been corrected.\n");
      allValuesInRange = false;
    }
    if (g.tau_s > maxTau_s)
    {
      g.tau_s = maxTau_s;
      printf("Gesture time constant was too high and has been corrected.\n");
      allValuesInRange = false;
    }

    // **************************************************************
    
    if (gestureNode->getAttributeInt("neutral") > 0)
    {
      g.neutral = true;
    }
    else
    {
      g.neutral = false;
    }

    appendGesture(g);
  }

  return true;
}


// ****************************************************************************
// ****************************************************************************

void GestureSequence::initGestureParams(Gesture &g)
{
  g.duration_s = 0.1;
  g.sVal = "";
  g.slope = 0.0;
  if (nominalValues)
  {
    g.dVal = 0.0;
  }
  else
  {
    g.dVal = 0.5*(minValue + maxValue);
  }

  g.neutral = false;
  g.tau_s = GesturalScore::DEFAULT_TIME_CONSTANT_S;
}


// ****************************************************************************
// ****************************************************************************

void GestureSequence::limitGestureParams(Gesture &g)
{
  // ****************************************************************

  if (g.duration_s < 0.001*MIN_GESTURE_DURATION_MS)
  {
    g.duration_s = 0.001*MIN_GESTURE_DURATION_MS;
  }

  if (g.duration_s > 0.001*MAX_GESTURE_DURATION_MS)
  {
    g.duration_s = 0.001*MAX_GESTURE_DURATION_MS;
  }

  // ****************************************************************

  if (g.dVal < minValue)
  {
    g.dVal = minValue;
  }

  if (g.dVal > maxValue)
  {
    g.dVal = maxValue;
  }

  // ****************************************************************

  if (g.slope < minSlope)
  {
    g.slope = minSlope;
  }

  if (g.slope > maxSlope)
  {
    g.slope = maxSlope;
  }

  // ****************************************************************

  if (g.tau_s < minTau_s)
  {
    g.tau_s = minTau_s;
  }

  if (g.tau_s > maxTau_s)
  {
    g.tau_s = maxTau_s;
  }
}


// ****************************************************************************
/// Constructor. Init the variables.
// ****************************************************************************

GesturalScore::GesturalScore(VocalTract *vocalTract, Glottis *glottis)
{
  this->vocalTract = vocalTract;
  this->glottis = glottis;

  pos = 0;
  leftTube = new Tube();
  rightTube = new Tube();
  leftTubeIndex = -1;


  // ****************************************************************
  // Init the resulting curve signal vectors.
  // ****************************************************************

  int i;
  
  // Reserve the capacity for the maximum number of samples.
  // The actual size will be adjusted depending on the length
  // of the gestural score.

  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    tractParamCurve[i].resize(MAX_CURVE_SAMPLES);
    tractParamTargets[i].resize(MAX_PARAM_TARGETS);
  }

  for (i=0; i < Glottis::MAX_CONTROL_PARAMS; i++)
  {
    glottisParamCurve[i].resize(MAX_CURVE_SAMPLES);
    glottisParamTargets[i].resize(MAX_PARAM_TARGETS);
  }

  inputSignal.resize(MAX_CURVE_SAMPLES);
  tauSignal.resize(MAX_CURVE_SAMPLES);

  // ****************************************************************
  // Init the gesture sequences.
  // ****************************************************************

  gestures[VOWEL_GESTURE].init("Vowel gestures", "vowel-gestures", "", 
    0.0, 0.0, 0.0, 0.0, 0.005, 0.04, true);

  gestures[LIP_GESTURE].init("Lip gestures", "lip-gestures", "", 
    0.0, 0.0, 0.0, 0.0, 0.005, 0.04, true);

  gestures[TONGUE_TIP_GESTURE].init("Tongue tip g.", "tongue-tip-gestures", "", 
    0.0, 0.0, 0.0, 0.0, 0.005, 0.04, true);

  gestures[TONGUE_BODY_GESTURE].init("Tongue body g.", "tongue-body-gestures", "", 
    0.0, 0.0, 0.0, 0.0, 0.005, 0.04, true);
  
  // The min and max values MUST correspond to the VO parameter in the vocal tract model!
  gestures[VELIC_GESTURE].init("Velic gestures", "velic-gestures", "", 
    -0.1, 1.0, 0.0, 0.0, 0.01, 0.04, false);

  gestures[GLOTTAL_SHAPE_GESTURE].init("Glottal shape g.", "glottal-shape-gestures", "", 
    0.0, 0.0, 0.0, 0.0, 0.01, 0.04, true);

  gestures[F0_GESTURE].init("F0 gestures", "f0-gestures", "st",
    48.0, 108.0, -80.0, 80.0, 0.005, 0.04, false);

  // The min and max values should correspond to the pressure parameter in the glottis models!
  gestures[PRESSURE_GESTURE].init("Lung pressure g.", "lung-pressure-gestures", "dPa", 
    0.0, 16000.0, 0.0, 0.0, 0.005, 0.04, false);

  initTestScore();
}


// ****************************************************************************
/// Destructor.
// ****************************************************************************

GesturalScore::~GesturalScore()
{
  // Do nothing...
}


// ****************************************************************************
/// Clear the gestural score.
// ****************************************************************************

void GesturalScore::clear()
{
  int i;

  pos = 0;
  leftTubeIndex = -1;

  for (i=0; i < NUM_GESTURE_TYPES; i++)
  {
    gestures[i].clear();
  }

  calcCurves();
}


// ****************************************************************************
/// Create a test gestural score.
// ****************************************************************************

void GesturalScore::initTestScore()
{
  Gesture g;
  int i;

  // Delete all potential previous gestures.
  for (i = 0; i < NUM_GESTURE_TYPES; i++)
  {
    gestures[i].clear();
  }

  // One vowel gesture
  
  g.duration_s = 0.6;
  g.dVal = 0.0;
  g.slope = 0.0;
  g.sVal = "a";
  g.neutral = false;
  g.tau_s = DEFAULT_TIME_CONSTANT_S;
  gestures[VOWEL_GESTURE].appendGesture(g);

  // F0-gestures
  
  g.duration_s = 0.01;
  g.dVal = 80.0;
  g.slope = 0.0;
  g.sVal = "";
  g.neutral = false;
  g.tau_s = 0.03;
  gestures[F0_GESTURE].appendGesture(g);
  
  g.duration_s = 0.3;
  g.dVal = 83.0;
  g.slope = 0.0;
  g.tau_s = 0.03;
  gestures[F0_GESTURE].appendGesture(g);

  g.duration_s = 0.3;
  g.dVal = 78.0;
  g.slope = 0.0;
  g.tau_s = 0.03;
  gestures[F0_GESTURE].appendGesture(g);

  // Lung pressure gestures

  g.duration_s = 0.01;
  g.dVal = 0.0;
  g.slope = 0.0;
  g.sVal = "";
  g.neutral = false;
  g.tau_s = 0.005;
  gestures[PRESSURE_GESTURE].appendGesture(g);

  g.duration_s = 0.5;
  g.dVal = 8000.0;
  g.tau_s = 0.005;
  gestures[PRESSURE_GESTURE].appendGesture(g);

  g.duration_s = 0.1;
  g.dVal = 0.0;
  g.tau_s = 0.005;
  gestures[PRESSURE_GESTURE].appendGesture(g);

  // Glottal shape gesture
  
  g.duration_s = 0.6;
  g.dVal = 0.0;
  g.slope = 0.0;
  g.sVal = "modal";
  g.neutral = false;
  g.tau_s = DEFAULT_TIME_CONSTANT_S;
  gestures[GLOTTAL_SHAPE_GESTURE].appendGesture(g);

}


// ****************************************************************************
/// Create the gestural score from the given segment sequence (i.e., phones
/// and their durations).
// ****************************************************************************

void GesturalScore::createFromSegmentSequence(SegmentSequence *origSegmentSequence)
{
  const double MIN_GESTURE_DURATION_S = (double)MIN_GESTURE_DURATION_MS * 0.001;
  int i;

  printf("Creating a gestural score from %d segments.\n", origSegmentSequence->numSegments());


  // ****************************************************************
  // Delete all potential previous gestures.
  // ****************************************************************

  gestures[VOWEL_GESTURE].clear();
  gestures[LIP_GESTURE].clear();
  gestures[TONGUE_TIP_GESTURE].clear();
  gestures[TONGUE_BODY_GESTURE].clear();
  gestures[VELIC_GESTURE].clear();
  gestures[GLOTTAL_SHAPE_GESTURE].clear();
  gestures[PRESSURE_GESTURE].clear();

  // Do NOT delete the F0 gestures, as they are not affected by the segment sequence.

  // ****************************************************************
  // Pre-process the segment sequence so that it is plain SAMPA:
  // * Unify all kinds of pause labels to the empty string.
  // * Split affricatives into the plosive and the fricative.
  // ****************************************************************

  int numSegments = origSegmentSequence->numSegments();
  Segment *origSegment;
  string n;     // The segment name
  double d;

  // The target segment sequence.
  SegmentSequence segmentSequence;

  for (i = 0; i < numSegments; i++)
  {
    origSegment = origSegmentSequence->getSegment(i);
    n = origSegment->value[Segment::NAME_INDEX];
    d = origSegment->duration_s;

    if (Sampa::isPhoneme(n))
    {
      // Split affricates into a plosive and a fricative.
      if ((n == "pf") || (n == "ts") || (n == "tS") || (n == "dZ"))
      {
        // The duration ratio of the plosive has been determined
        // from audio recordings of Peter Birkholz.
        double plosiveRatio = 0.4;
        if (n == "pf") { plosiveRatio = 0.4;  } 
        else
          if (n == "ts") { plosiveRatio = 0.33; }
          else
            if (n == "tS") { plosiveRatio = 0.33; }
            else
              if (n == "dZ") { plosiveRatio = 0.4; }

        segmentSequence.appendSegment(n.substr(0, 1), plosiveRatio * d);
        segmentSequence.appendSegment(n.substr(1, 1), (1.0 - plosiveRatio) * d);
      }
      else
      {
        segmentSequence.appendSegment(n, d);
      }
    }
    else
    {
      // Everything that is not a valid SAMPA phoneme is appended as
      // a pause (empty string).
      segmentSequence.appendSegment("", d);
    }
  }


  // ****************************************************************
  // Determine for the alveolars /d,t,n,l/ whether they should be
  // realized as *post*alveolar consonants. This is always the case
  // when they are next to /S,Z/ in the same consonant cluster.
  // ****************************************************************

  numSegments = segmentSequence.numSegments();
  Segment* s;

  // Initialize all entries as false.
  vector<bool> isPostalveolar(numSegments, false);

  for (i = 0; i < numSegments; i++)
  {
    s = segmentSequence.getSegment(i);
    n = s->value[Segment::NAME_INDEX];

    // The i-th phone is a postalveolar fricative.
    if ((n == "S") || (n == "Z"))
    {
      isPostalveolar[i] = true;

      // Look at maximal 3 phones back.
      if ((i - 1 >= 0) && (Sampa::isAlveolar(segmentSequence.getSegment(i - 1)->value[Segment::NAME_INDEX])))
      {
        isPostalveolar[i - 1] = true;

        if ((i - 2 >= 0) && (Sampa::isAlveolar(segmentSequence.getSegment(i - 2)->value[Segment::NAME_INDEX])))
        {
          isPostalveolar[i - 2] = true;

          if ((i - 3 >= 0) && (Sampa::isAlveolar(segmentSequence.getSegment(i - 3)->value[Segment::NAME_INDEX])))
          {
            isPostalveolar[i - 3] = true;
          }
        }
      }

      // Look at maximal 3 phones ahead.
      if ((i + 1 < numSegments) && (Sampa::isAlveolar(segmentSequence.getSegment(i + 1)->value[Segment::NAME_INDEX])))
      {
        isPostalveolar[i + 1] = true;

        if ((i + 2 < numSegments) && (Sampa::isAlveolar(segmentSequence.getSegment(i + 2)->value[Segment::NAME_INDEX])))
        {
          isPostalveolar[i + 2] = true;

          if ((i + 3 < numSegments) && (Sampa::isAlveolar(segmentSequence.getSegment(i + 3)->value[Segment::NAME_INDEX])))
          {
            isPostalveolar[i + 3] = true;
          }
        }
      }
    }
  }

  // ****************************************************************
  // Run through all segments and find the beginning and the end of 
  // the valid (non-pause) segments, respectively.
  // ****************************************************************

  numSegments = segmentSequence.numSegments();
  double segmentBegin_s = 0.0;
  double segmentEnd_s = 0.0;

  double firstValidSegmentBegin_s = 1000000.0;
  double lastValidSegmentEnd_s = 0.0;
  string firstValidSegmentName = "";
  string lastValidSegmentName = "";
  int firstValidSegmentIndex = 0;
  int lastValidSegmentIndex = 0;

  for (i = 0; i < numSegments; i++)
  {
    s = segmentSequence.getSegment(i);
    n = s->value[Segment::NAME_INDEX];
    segmentBegin_s = segmentSequence.getSegmentBegin_s(i);
    segmentEnd_s = segmentSequence.getSegmentEnd_s(i);

    // Skip empty segments (pauses).

    if (n.empty())
    {
      continue;
    }
    
    // All other segments are valid phonemes.

    if (segmentBegin_s < firstValidSegmentBegin_s)
    {
      firstValidSegmentBegin_s = segmentBegin_s;
      firstValidSegmentName = s->value[Segment::NAME_INDEX];
      firstValidSegmentIndex = i;
    }
    if (segmentEnd_s > lastValidSegmentEnd_s)
    {
      lastValidSegmentEnd_s = segmentEnd_s;
      lastValidSegmentName = s->value[Segment::NAME_INDEX];
      lastValidSegmentIndex = i;
    }
  }

  if (lastValidSegmentEnd_s < MIN_GESTURE_DURATION_S)
  {
    return;
  }


  // ****************************************************************
  // Create the basic lung gesture from the beginning of the first
  // valid segment to the end of the last valid segment.
  // ****************************************************************

  Gesture g;
  double duration_s = firstValidSegmentBegin_s;

  // If the first valid segment is a plosive, start the raise
  // of subglottal pressure 50 ms before the segment beginning.
  // For all other segments, start 20 ms before the beginning
  // of the first segment.

  if (Sampa::isPlosive(firstValidSegmentName))
  {
    duration_s -= 0.050;
  }
  else
  {
    duration_s -= 0.020;
  }

  if (duration_s < MIN_GESTURE_DURATION_S)
  {
    duration_s = MIN_GESTURE_DURATION_S;
  }

  // Zero pressure gesture at the beginning.
  g.duration_s = duration_s;
  g.dVal = 0.0;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = "";
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[PRESSURE_GESTURE].appendGesture(g);

  
  // If the final segment is a plosive, the lung pressure should 
  // start to drop 70 ms before the end of the segment in order to
  // get a good audible plosive release burst.
  // For all other sounds, the lung pressure should start to drop
  // already 120 ms before the end of the last segment, so that it 
  // is close to zero at the end of the last segment.

  double lungGestureEnd_s = lastValidSegmentEnd_s - 0.120;

  if ((lastValidSegmentName == "p") || 
    (lastValidSegmentName == "k"))
  {
    lungGestureEnd_s = lastValidSegmentEnd_s - 0.090;
  }
  else
  if (lastValidSegmentName == "t")
  {
    lungGestureEnd_s = lastValidSegmentEnd_s - 0.070;
  }

   
  duration_s = lungGestureEnd_s - gestures[PRESSURE_GESTURE].getDuration_s();
  if (duration_s > MIN_GESTURE_DURATION_S)
  {
    // The 800 Pa gesture
    g.duration_s = duration_s;
    g.dVal = 8000.0;      // = 800 Pa
    g.neutral = false;
    g.slope = 0.0;
    g.sVal = "";
    g.tau_s = 0.005;    // Only 5 ms

    gestures[PRESSURE_GESTURE].appendGesture(g);
  }

  // Pressure goes back to zero.
  g.duration_s = 0.17;    // 170 ms
  g.dVal = 0.0;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = "";
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[PRESSURE_GESTURE].appendGesture(g);

  // The lung pressure tier is now finished!

  // ****************************************************************
  // Create a gesture for modal phonation up to the end of the last
  // valid segment. Other glottal gestures for consonants may later
  // overwrite the modal gesture locally.
  // ****************************************************************

  g.duration_s = lastValidSegmentEnd_s;
  g.dVal = 0.0;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = "modal";
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[GLOTTAL_SHAPE_GESTURE].appendGesture(g);


  // ****************************************************************
  // Run through all segments and create the VOWEL gestures.
  // ****************************************************************

  bool isVowel = false;
  bool isDiphthong = false;
  string vowelName;         // Name of the vowel gesture
  string diphthongName[2];   // Names of the start and end gestures
  string segmentName;

  for (i = 0; i < numSegments; i++)
  {
    s = segmentSequence.getSegment(i);
    segmentName = s->value[Segment::NAME_INDEX];
    segmentBegin_s = segmentSequence.getSegmentBegin_s(i);
    segmentEnd_s = segmentSequence.getSegmentEnd_s(i);

    if ((Sampa::isVowel(segmentName) == false) && (Sampa::isDiphthong(segmentName) == false))
    {
      continue;
    }

    isVowel = false;
    isDiphthong = false;

    // Primary diphthongs must be searched for first.

    if (segmentName == "aI")
    {
      isDiphthong = true;
      diphthongName[0] = "a";
      diphthongName[1] = "e";
    }
    
    if (segmentName == "aU")
    {
      isDiphthong = true;
      diphthongName[0] = "a";
      diphthongName[1] = "o";
    }
    
    if (segmentName == "OY")
    {
      isDiphthong = true;
      diphthongName[0] = "O";
      diphthongName[1] = "e";
    }

    // Secondary diphthongs are searched for next.
    // The end vowels (one of two Tiefschwas) are selected according
    // to our study on secundary diphtongues (submitted to Journal of 
    // Phonetics, 2020)

    if ((segmentName == "i:6") || (segmentName == "i6"))
    {
      isDiphthong = true;
      diphthongName[0] = "i";
      diphthongName[1] = "6_mid";
    }
    
    if (segmentName == "I6")
    {
      isDiphthong = true;
      diphthongName[0] = "I";
      diphthongName[1] = "6_mid";
    }
    
    if ((segmentName == "y:6") || (segmentName == "y6"))
    {
      isDiphthong = true;
      diphthongName[0] = "y";
      diphthongName[1] = "6_mid";
    }
    
    if (segmentName == "Y6")
    {
      isDiphthong = true;
      diphthongName[0] = "Y";
      diphthongName[1] = "6_mid";
    }
    
    if ((segmentName == "e:6") || (segmentName == "e6"))
    {
      isDiphthong = true;
      diphthongName[0] = "e";
      diphthongName[1] = "6_low";
    }

    if (segmentName == "E6")
    {
      isDiphthong = true;
      diphthongName[0] = "E";
      diphthongName[1] = "6_low";
    }

    if (segmentName == "E:6")
    {
      isDiphthong = true;
      diphthongName[0] = "E:";
      diphthongName[1] = "6_low";
    }

    if ((segmentName == "2:6") || (segmentName == "26"))
    {
      isDiphthong = true;
      diphthongName[0] = "2";
      diphthongName[1] = "6_low";
    }
 
    if (segmentName == "96")
    {
      isDiphthong = true;
      diphthongName[0] = "9";
      diphthongName[1] = "6_low";
    }

    if (segmentName == "a:6")
    {
      isDiphthong = true;
      diphthongName[0] = "a";
      diphthongName[1] = "6_low";
    }

    if (segmentName == "a6")
    {
      isDiphthong = true;
      diphthongName[0] = "a";
      diphthongName[1] = "6_low";
    }

    if ((segmentName == "u:6") || (segmentName == "u6"))
    {
      isDiphthong = true;
      diphthongName[0] = "u";
      diphthongName[1] = "6_mid";
    }

    if (segmentName == "U6")
    {
      isDiphthong = true;
      diphthongName[0] = "U";
      diphthongName[1] = "6_mid";
    }

    if ((segmentName == "o:6") || (segmentName == "o6"))
    {
      isDiphthong = true;
      diphthongName[0] = "o";
      diphthongName[1] = "6_low";
    }

    if (segmentName == "O6")
    {
      isDiphthong = true;
      diphthongName[0] = "O";
      diphthongName[1] = "6_low";
    }

    // Now we search for all other vowels.

    if ((segmentName == "a:") || (segmentName == "a"))
    {
      isVowel = true;
      vowelName = "a";
    }

    if ((segmentName == "e:") || (segmentName == "e"))
    {
      isVowel = true;
      vowelName = "e";
    }

    if ((segmentName == "i:") || (segmentName == "i"))
    {
      isVowel = true;
      vowelName = "i";
    }

    if (segmentName == "I")
    {
      isVowel = true;
      vowelName = "I";
    }

    if ((segmentName == "o:") || (segmentName == "o"))
    {
      isVowel = true;
      vowelName = "o";
    }

    if (segmentName == "O")
    {
      isVowel = true;
      vowelName = "O";
    }

    if ((segmentName == "u:") || (segmentName == "u"))
    {
      isVowel = true;
      vowelName = "u";
    }

    if (segmentName == "U")
    {
      isVowel = true;
      vowelName = "U";
    }

    if (segmentName == "E:")
    {
      isVowel = true;
      vowelName = "E:";
    }

    if (segmentName == "E")
    {
      isVowel = true;
      vowelName = "E";
    }

    if ((segmentName == "2:") || (segmentName == "2"))
    {
      isVowel = true;
      vowelName = "2";
    }

    if (segmentName == "9")
    {
      isVowel = true;
      vowelName = "9";
    }

    if ((segmentName == "y:") || (segmentName == "y"))
    {
      isVowel = true;
      vowelName = "y";
    }

    if (segmentName == "Y")
    {
      isVowel = true;
      vowelName = "Y";
    }

    if (segmentName == "@")
    {
      isVowel = true;
      vowelName = "@";
    }

    if (segmentName == "6")
    {
      isVowel = true;
      vowelName = "6_low";
    }

    // **************************************************************
    // Append a vowel gesture to the vowel tier.
    // **************************************************************

    if (isVowel)
    {
      double tierDuration_s = gestures[VOWEL_GESTURE].getDuration_s();
      double gestureEnd_s = segmentEnd_s - 0.020;
      double gestureDuration_s = gestureEnd_s - tierDuration_s;
      if (gestureDuration_s > MIN_GESTURE_DURATION_S)
      {
        g.duration_s = gestureDuration_s;
        g.dVal = 0.0;
        g.neutral = false;
        g.slope = 0.0;
        g.sVal = vowelName;
        g.tau_s = DEFAULT_TIME_CONSTANT_S;

        gestures[VOWEL_GESTURE].appendGesture(g);
      }
    }
    else

    // **************************************************************
    // Append two diphthong gestures to the vowel tier.
    // **************************************************************

    if (isDiphthong)
    {
      double tierDuration_s = gestures[VOWEL_GESTURE].getDuration_s();
      
      // The start of the 2nd vowel gesture seems to depend on the 
      // duration of the diphthong. For short diphthongs (< 130 ms),
      // the 2nd vowel gesture starts right at the beginning of the
      // diphthong. For longer diphthongs, the 2nd vowel gesture 
      // should start up to 70 ms later (for both primary and 
      // secondary diphthongs):
      // These values were determined subjectively so far!
      // This needs to be formally evaluated in a perception experiment.

      const double SHORT_DIPHTHONG_DURATION_S = 0.130;    // 130 ms
      const double MAX_LAG_S = 0.070;   // 70 ms
      double lag_s = 0.0;

      double segmentDuration_s = segmentEnd_s - segmentBegin_s;
      if (segmentDuration_s > SHORT_DIPHTHONG_DURATION_S)
      {
        lag_s = segmentDuration_s - SHORT_DIPHTHONG_DURATION_S;
        if (lag_s > MAX_LAG_S)
        {
          lag_s = MAX_LAG_S;
        }
      }

      double firstGestureEnd_s = segmentBegin_s + lag_s;
      double secondGestureEnd_s = segmentEnd_s - 0.020;

      // Add the first diphthong gesture.

      double gestureDuration_s = firstGestureEnd_s - tierDuration_s;
      if (gestureDuration_s > MIN_GESTURE_DURATION_S)
      {
        g.duration_s = gestureDuration_s;
        g.dVal = 0.0;
        g.neutral = false;
        g.slope = 0.0;
        g.sVal = diphthongName[0];
        g.tau_s = DEFAULT_TIME_CONSTANT_S;

        gestures[VOWEL_GESTURE].appendGesture(g);
      }
      
      // Add the second diphthong gesture.

      tierDuration_s = gestures[VOWEL_GESTURE].getDuration_s();
      gestureDuration_s = secondGestureEnd_s - tierDuration_s;
      if (gestureDuration_s > MIN_GESTURE_DURATION_S)
      {
        g.duration_s = gestureDuration_s;
        g.dVal = 0.0;
        g.neutral = false;
        g.slope = 0.0;
        g.sVal = diphthongName[1];
        g.tau_s = DEFAULT_TIME_CONSTANT_S;

        gestures[VOWEL_GESTURE].appendGesture(g);
      }
    }

  }     // End of segment loop.


  // ****************************************************************
  // Add a final vowel schwa gesture that fills the gap until the end
  // of the segment sequence.
  // ****************************************************************

  double tierDuration_s = gestures[VOWEL_GESTURE].getDuration_s();
  double segmentSequenceDuration_s = segmentSequence.getDuration_s();
  double gestureDuration_s = segmentSequenceDuration_s - tierDuration_s;

  if (gestureDuration_s > MIN_GESTURE_DURATION_S)
  {
    g.duration_s = gestureDuration_s;
    g.dVal = 0.0;
    g.neutral = false;
    g.slope = 0.0;
    g.sVal = "@";
    g.tau_s = DEFAULT_TIME_CONSTANT_S;

    gestures[VOWEL_GESTURE].appendGesture(g);
  }

  // The vowel tier is now finished!


  // ****************************************************************
  // Run through all segments and create the FRICATIVE gestures.
  // ****************************************************************

  for (i = 0; i < numSegments; i++)
  {
    s = segmentSequence.getSegment(i);
    segmentName = s->value[Segment::NAME_INDEX];
    segmentBegin_s = segmentSequence.getSegmentBegin_s(i);
    segmentEnd_s = segmentSequence.getSegmentEnd_s(i);

    bool isFricative = false;
    string oralGestureName;
    string glottalGestureName;
    GestureType gestureType = LIP_GESTURE;

    // **************************************************************
    // Determine the segment.
    // **************************************************************

    // Search for all FRICATIVES.

    if (segmentName == "f")
    {
      isFricative = true;
      gestureType = LIP_GESTURE;
      oralGestureName = "ll-dental-fricative";
      glottalGestureName = "voiceless-fricative";
    }

    if (segmentName == "v")
    {
      isFricative = true;
      gestureType = LIP_GESTURE;
      oralGestureName = "ll-dental-fricative";
      glottalGestureName = "voiced-fricative";
    }

    if (segmentName == "T")
    {
      isFricative = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-dental-fricative";
      glottalGestureName = "voiceless-fricative";
    }

    if (segmentName == "D")
    {
      isFricative = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-dental-fricative";
      glottalGestureName = "voiced-fricative";
    }

    if (segmentName == "s")
    {
      isFricative = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-alveolar-fricative";
      glottalGestureName = "voiceless-fricative";
    }

    if (segmentName == "z")
    {
      isFricative = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-alveolar-fricative";
      glottalGestureName = "voiced-fricative";
    }

    if (segmentName == "S")
    {
      isFricative = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-postalveolar-fricative";
      glottalGestureName = "voiceless-fricative";
    }

    if (segmentName == "Z")
    {
      isFricative = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-postalveolar-fricative";
      glottalGestureName = "voiced-fricative";
    }

    if (segmentName == "C")
    {
      isFricative = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-palatal-fricative";
      glottalGestureName = "voiceless-fricative";
    }

    if (segmentName == "j")
    {
      isFricative = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-palatal-fricative";
      // Take modal phonation here, otherwise the noise gets too strong!
      glottalGestureName = "modal";
    }

    if (segmentName == "x")
    {
      isFricative = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-uvular-fricative";
      glottalGestureName = "voiceless-fricative";
    }

    if ((segmentName == "r") || (segmentName == "R"))
    {
      isFricative = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-uvular-fricative";
      glottalGestureName = "voiced-fricative";
    }

    // **************************************************************
    // Add the gestures for the fricative.
    // **************************************************************

    if (isFricative)
    {
      // This delay was precisely tuned so that the degree of glottal
      // abduction is the same at the segment beginning and end.

      const double DELTA_S = 0.057;

      // Add the oral closing gesture.
      // Default values for the voiceless fricatives:

      double oralGestureBegin_s = segmentBegin_s - DELTA_S - 0.015;
      double oralGestureEnd_s = segmentEnd_s - DELTA_S + 0.010;

      // Is it a voiced fricative?
      if (Sampa::isVoiced(segmentName))
      {
        oralGestureBegin_s = segmentBegin_s - DELTA_S - 0.020;
      }

      // If this fricative is the very first segment, start the
      // oral gesture at the very beginning.

      if (i == firstValidSegmentIndex)
      {
        oralGestureBegin_s = 0.0;
      }

      g.duration_s = oralGestureEnd_s - oralGestureBegin_s;
      g.dVal = 0.0;
      g.neutral = false;
      g.slope = 0.0;
      g.sVal = oralGestureName;
      g.tau_s = DEFAULT_TIME_CONSTANT_S;

      gestures[gestureType].putGesture(g, oralGestureBegin_s);

      // Add the corresponding glottal gesture. It must have the same
      // length as the segment, so that a cluster of fricatives has
      // a continuously open glottis. But it starts earlier than the 
      // segment, so that the glottal opening is the same at the 
      // beginning and end of the segment.

      double glottalGestureBegin_s = segmentBegin_s - DELTA_S;
      double glottalGestureEnd_s = segmentEnd_s - DELTA_S;

      // If this fricative is the very first segment, start the
      // glottal gesture at the very beginning.

      if (i == firstValidSegmentIndex)
      {
        glottalGestureBegin_s = 0.0;
      }

      // If the current segment is a VOICELESS fricative and
      // the next segment is a VOICELESS plosive, extend the 
      // glottal gesture by 30 ms (empirical value).

      Segment *nextSegment = segmentSequence.getSegment(i + 1);
      string nextSegmentName = "";
      if (nextSegment != NULL)
      {
        nextSegmentName = nextSegment->value[Segment::NAME_INDEX];
      }

      if ((Sampa::isVoiced(segmentName) == false) &&
        ((nextSegmentName == "p") || (nextSegmentName == "t") || (nextSegmentName == "k")))
      {
        glottalGestureEnd_s += 0.030;
      }
                     
      if (glottalGestureEnd_s >= glottalGestureBegin_s + MIN_GESTURE_DURATION_S)
      {
        g.duration_s = glottalGestureEnd_s - glottalGestureBegin_s;
        g.dVal = 0.0;
        g.neutral = false;
        g.slope = 0.0;
        g.sVal = glottalGestureName;
        g.tau_s = DEFAULT_TIME_CONSTANT_S;

        gestures[GLOTTAL_SHAPE_GESTURE].putGesture(g, glottalGestureBegin_s);
      }
    }     // End of (isFricative == true)

  }   // End of segment loop

  // The fricatives are now finished.


  // ****************************************************************
  // Run through all segments and create the PLOSIVES, NASALS, the
  // LATERAL, the GLOTTAL STOP, and the GLOTTAL FRICATIVE.
  // It is important to create the gestures with a precise 
  // oral closure after the fricatives and the vowels.
  // ****************************************************************

  string prevOralGestureName = "";

  for (i = 0; i < numSegments; i++)
  {
    s = segmentSequence.getSegment(i);
    segmentName = s->value[Segment::NAME_INDEX];
    segmentBegin_s = segmentSequence.getSegmentBegin_s(i);
    segmentEnd_s = segmentSequence.getSegmentEnd_s(i);

    string prevSegmentName = "";
    string nextSegmentName = "";
    
    if (segmentSequence.getSegment(i - 1) != NULL)
    {
      prevSegmentName = segmentSequence.getSegment(i - 1)->value[Segment::NAME_INDEX];
    }
    
    if (segmentSequence.getSegment(i + 1) != NULL)
    {
      nextSegmentName = segmentSequence.getSegment(i + 1)->value[Segment::NAME_INDEX];
    }

    bool isPlosive = false;
    bool isNasal = false;
    string oralGestureName;
    string glottalGestureName = "";
    double voiceOnsetTime_s = 0.0;
    GestureType gestureType = LIP_GESTURE;

    // **************************************************************
    // Determine the current segment.
    // **************************************************************

    // Search for all PLOSIVES.
    // The voice onset times were taken from Klatt (1975):
    // "Voice onset time, frication, and aspiration ..."

    if (segmentName == "p")
    {
      isPlosive = true;
      gestureType = LIP_GESTURE;
      oralGestureName = "ll-labial-closure";
      glottalGestureName = "voiceless-plosive";
      if (Sampa::isFricative(nextSegmentName))
      {
        voiceOnsetTime_s = 0.010;
      }
      else
      {
        voiceOnsetTime_s = 0.047;
      }
    }

    if (segmentName == "b")
    {
      isPlosive = true;
      gestureType = LIP_GESTURE;
      oralGestureName = "ll-labial-closure";
      glottalGestureName = "voiced-plosive";
      voiceOnsetTime_s = 0.011;
    }

    if (segmentName == "t")
    {
      isPlosive = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-alveolar-closure";
      if (isPostalveolar[i])    // there is a /S,Z/ in the neighborhood
      {
        oralGestureName = "tt-postalveolar-closure";
      }
      glottalGestureName = "voiceless-plosive";
      if (Sampa::isFricative(nextSegmentName))
      {
        voiceOnsetTime_s = 0.015;
      }
      else
      {
        voiceOnsetTime_s = 0.050;
      }
    }

    if (segmentName == "d")
    {
      isPlosive = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-alveolar-closure";
      if (isPostalveolar[i])    // there is a /S,Z/ in the neighborhood
      {
        oralGestureName = "tt-postalveolar-closure";
      }
      glottalGestureName = "voiced-plosive";
      voiceOnsetTime_s = 0.017;
    }

    if (segmentName == "k")
    {
      isPlosive = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-velar-closure";
      glottalGestureName = "voiceless-plosive";
      if (Sampa::isFricative(nextSegmentName))
      {
        voiceOnsetTime_s = 0.020;
      }
      else
      {
        voiceOnsetTime_s = 0.070;
      }
    }

    if (segmentName == "g")
    {
      isPlosive = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-velar-closure";
      glottalGestureName = "voiced-plosive";
      voiceOnsetTime_s = 0.027;
    }

    // Check the nasals.

    if (segmentName == "m")
    {
      isNasal = true;
      gestureType = LIP_GESTURE;
      oralGestureName = "ll-labial-closure";
    }

    if (segmentName == "n")
    {
      isNasal = true;
      gestureType = TONGUE_TIP_GESTURE;
      oralGestureName = "tt-alveolar-closure";
      if (isPostalveolar[i])    // there is a /S,Z/ in the neighborhood
      {
        oralGestureName = "tt-postalveolar-closure";
      }
    }

    if (segmentName == "N")
    {
      isNasal = true;
      gestureType = TONGUE_BODY_GESTURE;
      oralGestureName = "tb-velar-closure";
    }

    // **************************************************************
    // Add the gestures for the PLOSIVE.
    // **************************************************************

    if (isPlosive)
    {
      // Beginning and end of the closure:

      double closureBegin_s = segmentBegin_s;
      double closureEnd_s = segmentEnd_s - voiceOnsetTime_s;

      // If this plosive is the first valid segment, start the 
      // oral gesture at t = 0.

      if (i == firstValidSegmentIndex)
      {
        closureBegin_s = 0.0;
      }
      else

      // If this plosive is NOT the very first segment, make sure
      // that the closure duration is not be shorter than half the
      // segment duration and not shorter than 20 ms.
      {
        const double MIN_CLOSURE_DURATION_S = 0.020;
        double segmentDur_s = segmentEnd_s - segmentBegin_s;

        if (closureEnd_s - closureBegin_s < 0.5 * segmentDur_s)
        {
          closureEnd_s = closureBegin_s + 0.5 * segmentDur_s;
        }

        if (closureEnd_s - closureBegin_s < MIN_CLOSURE_DURATION_S)
        {
          closureEnd_s = closureBegin_s + MIN_CLOSURE_DURATION_S;
        }
      }

      // Beginning and end of the glottal gesture:

      double glottalGestureBegin_s = closureBegin_s - 0.010;    // 10 ms before the implosion
      double glottalGestureEnd_s = segmentEnd_s - 0.040;

      // If this plosive is the first valid segment, let the glottal
      // gesture start at t = 0.

      if (i == firstValidSegmentIndex)
      {
        glottalGestureBegin_s = 0.0;
      }

      // Shorten the glottal gesture if the next segment is a voiced fricative.
      if ((Sampa::isFricative(nextSegmentName)) && (Sampa::isVoiced(nextSegmentName)))
      {
        glottalGestureEnd_s -= 0.020;
      }

      // In the case that the plosive is followed by a nasal or 
      // lateral with the same place of articulation, the release
      // is right at the end of the plosive segment.

      if ((((segmentName == "p") || (segmentName == "b")) && (nextSegmentName == "m")) ||
        (((segmentName == "t") || (segmentName == "d")) && ((nextSegmentName == "n") || (nextSegmentName == "l"))) ||
        (((segmentName == "k") || (segmentName == "g")) && (nextSegmentName == "N")))
      {
        closureEnd_s = segmentEnd_s;
        // Make the glottal gesture much shorter.
        glottalGestureEnd_s = segmentEnd_s - 0.080;
      }

      // If the current segment is /p,b,k,g/ and the next segment
      // is /l/ (different place of articulation), shorten the 
      // glottal gesture by 15 ms.

      if ((nextSegmentName == "l") &&
        ((segmentName == "k") || (segmentName == "g") ||
        (segmentName == "p") || (segmentName == "b")))
      {
        glottalGestureEnd_s = segmentEnd_s - 0.055;
      }

      // Add the oral closing gesture.

      bool connectToPrevGesture = (oralGestureName == prevOralGestureName);
      addClosingGesture(gestureType, oralGestureName, closureBegin_s, closureEnd_s, connectToPrevGesture);

      // Add the corresponding glottal gesture (the velic port needs to be closed).

      // If this (voiceless) plosive is the last segment
      // extend the glottal gesture by 100 ms for a good release burst.
      
      if ((i == lastValidSegmentIndex) && (Sampa::isVoiced(segmentName) == false))
      {
        glottalGestureEnd_s += 0.100;
      }

      if (glottalGestureEnd_s - glottalGestureBegin_s < MIN_GESTURE_DURATION_S)
      {
        glottalGestureEnd_s = glottalGestureBegin_s + MIN_GESTURE_DURATION_S;
      }

      g.duration_s = glottalGestureEnd_s - glottalGestureBegin_s;
      g.dVal = 0.0;
      g.neutral = false;
      g.slope = 0.0;
      g.sVal = glottalGestureName;
      g.tau_s = DEFAULT_TIME_CONSTANT_S;

      gestures[GLOTTAL_SHAPE_GESTURE].putGesture(g, glottalGestureBegin_s);
    }
    else

    // **************************************************************
    // Add the gestures for the NASAL.
    // **************************************************************

    if (isNasal)
    {
      double closureBegin_s = segmentBegin_s;
      double closureEnd_s = segmentEnd_s;

      // Default for vowels.
      double velicOpeningBegin_s = segmentBegin_s - 0.030;
      double velicOpeningEnd_s = segmentEnd_s + 0.030;

      // Some special cases when plosives or fricatives are the
      // neighbouring segments.

      if (Sampa::isPlosive(prevSegmentName))
      {
        velicOpeningBegin_s = segmentBegin_s - 0.020;
      }
      if (Sampa::isPlosive(nextSegmentName))
      {
        velicOpeningEnd_s = segmentEnd_s;
      }

      if (Sampa::isFricative(prevSegmentName))
      {
        velicOpeningBegin_s = segmentBegin_s - 0.005;
      }
      if (Sampa::isFricative(nextSegmentName))
      {
        velicOpeningEnd_s = segmentEnd_s + 0.005;
      }

      if (prevSegmentName == "l")
      {
        velicOpeningBegin_s = segmentBegin_s - 0.035;
      }
      if (nextSegmentName == "l")
      {
        velicOpeningEnd_s = segmentEnd_s + 0.005;
      }

      // Add the oral closing gesture.

      bool connectToPrevGesture = (oralGestureName == prevOralGestureName);

      addClosingGesture(gestureType, oralGestureName, closureBegin_s, closureEnd_s, connectToPrevGesture);

      // Add the velic opening gesture.
      addVelicOpeningGesture(velicOpeningBegin_s, velicOpeningEnd_s);
    }

    // **************************************************************
    // Create the gestures for /h/, the glottal stop, and the lateral.
    // **************************************************************

    if (segmentName == "h")
    {
      double gestureBegin_s = segmentBegin_s - 0.060;
      double gestureEnd_s = segmentEnd_s - 0.060;

      g.duration_s = gestureEnd_s - gestureBegin_s;
      g.dVal = 0.0;
      g.neutral = false;
      g.slope = 0.0;
      g.sVal = "h";
      g.tau_s = DEFAULT_TIME_CONSTANT_S;

      gestures[GLOTTAL_SHAPE_GESTURE].putGesture(g, gestureBegin_s);
    }

    if (segmentName == "?")
    {
      double gestureBegin_s = segmentBegin_s - 0.060;
      double gestureEnd_s = segmentEnd_s - 0.060;
      
      // If this is the first valid segment, start the stop gesture
      // at t = 0.
      if (i == firstValidSegmentIndex)
      {
        gestureBegin_s = 0.0;
      }

      g.duration_s = gestureEnd_s - gestureBegin_s;
      g.dVal = 0.0;
      g.neutral = false;
      g.slope = 0.0;
      g.sVal = "stop";
      g.tau_s = DEFAULT_TIME_CONSTANT_S;

      gestures[GLOTTAL_SHAPE_GESTURE].putGesture(g, gestureBegin_s);
    }

    if (segmentName == "l")
    {
      // Make the gesture a little longer than the segment
      // and shift the gesture towards the left.
      double gestureBegin_s = segmentBegin_s - 0.055 - 0.010;
      double gestureEnd_s = segmentEnd_s - 0.055;

      g.duration_s = gestureEnd_s - gestureBegin_s;
      g.dVal = 0.0;
      g.neutral = false;
      g.slope = 0.0;
      g.sVal = "tt-alveolar-lateral";
      if (isPostalveolar[i])    // there is a /S,Z/ in the neighborhood
      {
        g.sVal = "tt-postalveolar-lateral";
      }
      g.tau_s = DEFAULT_TIME_CONSTANT_S;

      gestures[TONGUE_TIP_GESTURE].putGesture(g, gestureBegin_s);
    }

    // Keep in mind the last oral gesture name.
    prevOralGestureName = oralGestureName;

  }   // End of segment loop

}


// ****************************************************************************
/// Put a new closing gesture into the score and adjust its boundaries such
/// that the actual closure starts and ends at the given times.
/// When connectToPrevGesture = true, the beginning of the new gesture is set
/// to the end of the last non-neutral gesture.
// ****************************************************************************

void GesturalScore::addClosingGesture(GestureType gestureType, string gestureName, 
  double closureBegin_s, double closureEnd_s, bool connectToPrevGesture)
{
  const double TIME_STEP_S = 1.0 / (double)CURVE_SAMPLING_RATE;
  const double MAX_GESTURE_DURATION_S = 0.4;

  // Range check: The closure duration must be at least 20 ms.
  if (closureEnd_s < closureBegin_s + 0.020)
  {
    closureEnd_s = closureBegin_s + 0.020;
  }

  // Initial guesses for the gesture begin and end.
  double gestureBegin_s = closureBegin_s - 0.010;
  double gestureEnd_s = closureBegin_s;

  if (connectToPrevGesture)
  {
    // Make the start of the new gesture to the end of the last
    // previous non-neutral gesture (the type of gesture may differ).

    GestureSequence *sequence = &gestures[gestureType];
    int numGestures = sequence->numGestures();
    int prevValidGestureIndex = -1;
    int i;

    for (i = 0; i < numGestures; i++)
    {
      if ((sequence->getGesture(i)->neutral == false) &&
        (sequence->getGestureEnd_s(i) < closureBegin_s))
      {
        prevValidGestureIndex = i;
      }
    }

    if (prevValidGestureIndex != -1)
    {
      gestureBegin_s = sequence->getGestureEnd_s(prevValidGestureIndex);
    }
  }


  Gesture g;
  bool hasClosure = false;
  bool hasError = false;

  // The gesture must not begin at t < 0. 
  if (gestureBegin_s <= 0.0)
  {
    gestureBegin_s = 0.0;
    hasClosure = true;
  }

  // ****************************************************************
  // If connectToPrevGesture == false,
  // move the beginning of the gesture in small steps towards the left
  // until a closure is produced at the intended closure beginning.
  // ****************************************************************

  while ((connectToPrevGesture == false) && (hasClosure == false) && (hasError == false))
  {
    // Move the begin of the closing gesture a small step to the left.
    gestureBegin_s -= TIME_STEP_S;

    hasClosure = hasVocalTactClosure(gestureType, gestureName,
      gestureBegin_s, gestureEnd_s, closureBegin_s);

    if ((gestureBegin_s < 0.010) || (gestureEnd_s - gestureBegin_s > MAX_GESTURE_DURATION_S))
    {
      if (gestureEnd_s - gestureBegin_s > MAX_GESTURE_DURATION_S)
      {
        // Take the initial time for the gesture begin.
        gestureBegin_s = closureBegin_s - 0.010;
      }
      hasError = true;
    }
  }

  // ****************************************************************
  // If the vocal tract is closed at the intended end of the closure,
  // move the end of the gesture left in small steps, until the 
  // closure disappears. Otherwise,  move the end of the gesture 
  // right in small steps, until the closure appears.
  // ****************************************************************

  hasClosure = hasVocalTactClosure(gestureType, gestureName,
    gestureBegin_s, gestureEnd_s, closureEnd_s);

  if (hasClosure)
  {
    // Currently, the vocal tract is closed at the intended beginning
    // of the closure. We want it to stay like that.
    bool hasOnsetClosure = true;
    // At the moment, the vocal tract is closed at the intended end of 
    // the closure. We want it to open by moving the gesture boundary.
    bool hasOffsetClosure = true;
    bool gestureTooShort = false;

    while ((hasOnsetClosure) && (hasOffsetClosure) && (gestureTooShort == false))
    {
      // Move the end of the gesture a small step to the left.
      gestureEnd_s -= TIME_STEP_S;

      hasOnsetClosure = hasVocalTactClosure(gestureType, gestureName,
        gestureBegin_s, gestureEnd_s, closureBegin_s);

      hasOffsetClosure = hasVocalTactClosure(gestureType, gestureName,
        gestureBegin_s, gestureEnd_s, closureEnd_s);

      if (gestureEnd_s - gestureBegin_s < 0.010)
      {
        gestureTooShort = true;
      }
    }

    // If the gesture ending was moved too far towards the left (so 
    // that the vocal tract re-opened at the intended implosion time)
    // of if the gesture has become too short, undo the last step left.

    if ((hasOnsetClosure == false) || (gestureTooShort))
    {
      gestureEnd_s += TIME_STEP_S;
    }
  }
  else
  
  // ****************************************************************

  {
    hasError = false;
    while ((hasClosure == false) && (hasError == false))
    {
      // Move the end of the gesture a small step to the right.
      gestureEnd_s += TIME_STEP_S;

      hasClosure = hasVocalTactClosure(gestureType, gestureName,
        gestureBegin_s, gestureEnd_s, closureEnd_s);

      if (gestureEnd_s - gestureBegin_s > MAX_GESTURE_DURATION_S)
      {
        gestureEnd_s = gestureBegin_s + 0.010;
        hasError = true;
      }
    }
  }

  // ****************************************************************
  // Put the new gesture into the tier.
  // ****************************************************************

  g.duration_s = gestureEnd_s - gestureBegin_s;
  g.dVal = 0.0;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = gestureName;
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[gestureType].putGesture(g, gestureBegin_s);
}


// ****************************************************************************
/// This function returns true if the following condition is fulfilled
/// at testTime_s, when the given oral gesture is put into the score:
/// There is an oral closure with the primary articulator associated with
/// the gestureType (lips, tongue tip, or tongue back).
// ****************************************************************************

bool GesturalScore::hasVocalTactClosure(GestureType gestureType, string gestureName,
  double gestureBegin_s, double gestureEnd_s, double testTime_s)
{
  static GestureSequence storedGestures;
  static Tube tube;
  static Gesture g;
  static double tractParams[VocalTract::NUM_PARAMS];
  static double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  Tube::Section *ts = NULL;
  int i;
  bool hasClosure = false;

  // Store the old state of the gesture tier.

  storedGestures = gestures[gestureType];

  // ****************************************************************
  // Put the gesture into the tier.
  // ****************************************************************

  g.duration_s = gestureEnd_s - gestureBegin_s;
  g.dVal = 0.0;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = gestureName;
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[gestureType].putGesture(g, gestureBegin_s);

  // ****************************************************************
  // Calculate the parameter curves and get the vocal tract tube
  // state at the time of the intended closure.
  // ****************************************************************

  calcCurves();
  getParams(testTime_s, tractParams, glottisParams);

  vocalTract->setParams(tractParams);
  vocalTract->calculateAll();
  vocalTract->getTube(&tube);

  // ****************************************************************
  // Determine the position of the tongue tip (most anterior point).
  // ****************************************************************

  double tongueTipPos_cm = 0.0;
  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &tube.pharynxMouthSection[i];
    if (ts->articulator == Tube::TONGUE)
    {
      tongueTipPos_cm = ts->pos_cm + ts->length_cm;
    }
  }

  // ****************************************************************
  // Check whether the tube has a closure formed with the primary 
  // articulator given by the gesture type.
  // ****************************************************************

  // The distance from the tongue tip where the region of the
  // tongue back starts
  const double TONGUE_TIP_REGION_CM = 2.2;
  double pos_cm;

  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &tube.pharynxMouthSection[i];

    if (ts->area_cm2 < 2.0 * Tube::MIN_AREA_CM2)
    {
      // Position in the middle of the closed tube section.
      pos_cm = ts->pos_cm + 0.5 * ts->length_cm;

      if ((gestureType == LIP_GESTURE) && (ts->articulator == Tube::LOWER_LIP))
      {
        hasClosure = true;
      }

      if ((gestureType == TONGUE_TIP_GESTURE) && (ts->articulator == Tube::TONGUE) &&
        (pos_cm > tongueTipPos_cm - TONGUE_TIP_REGION_CM))
      {
        hasClosure = true;
      }

      if ((gestureType == TONGUE_BODY_GESTURE) && (ts->articulator == Tube::TONGUE) &&
        (pos_cm < tongueTipPos_cm - TONGUE_TIP_REGION_CM))
      {
        hasClosure = true;
      }
    }
  }

  // Restore the old state of the gesture tier.
  gestures[gestureType] = storedGestures;

  return hasClosure;
}


// ****************************************************************************
/// Adds a velic opening gesture such that the opening starts at openingBegin_s
// and ends at openingEnd_s.
// ****************************************************************************

void GesturalScore::addVelicOpeningGesture(double openingBegin_s, double openingEnd_s)
{
  const double TIME_STEP_S = 1.0 / (double)CURVE_SAMPLING_RATE;
  const double MAX_GESTURE_DURATION_S = 0.4;

  // Range check: The velic opening duration must be at least 20 ms.
  if (openingEnd_s < openingBegin_s + 0.020)
  {
    openingEnd_s = openingBegin_s + 0.020;
  }

  // ****************************************************************
  // Conservative initialization for the beginning and end of the 
  // gesture. 
  // For a target VO of 100%, this initialization gives  exactly 
  // the right gesture boundaries (for tau = 12 ms). 
  // For smaller values of the target VO, the gesture length will be
  // expanded to both sides below.
  // If the conext vowel is sightly nasalized, the initialized 
  // bondaries are not shifted anymore, because the velic port is 
  // already open at the intended beginning and end of the opening.
  // ****************************************************************

  double gestureBegin_s = openingBegin_s - 0.030;
  double gestureEnd_s = openingEnd_s - 0.090;

  // Range check: The gesture must not begin at t < 0. 
  if (gestureBegin_s <= 0.0)
  {
    gestureBegin_s = 0.0;
  }

  // Range check: The gesture end must be AFTER the beginning.
  if (gestureEnd_s < gestureBegin_s + 0.001)
  {
    gestureEnd_s = gestureBegin_s + 0.001;
  }

  // ****************************************************************
  // Move the beginning of the gesture in small steps towards the left
  // until an opening is produced at the intended opening beginning.
  // ****************************************************************

  Gesture g;
  bool hasOpening = false;
  bool hasError = false;

  hasOpening = hasVelicOpening(gestureBegin_s, gestureEnd_s, openingBegin_s);

  while ((hasOpening == false) && (hasError == false))
  {
    // Move the beginning of the closing gesture 5 ms to the left.
    gestureBegin_s -= TIME_STEP_S;

    hasOpening = hasVelicOpening(gestureBegin_s, gestureEnd_s, openingBegin_s);

    if ((gestureBegin_s < 0.010) || (gestureEnd_s - gestureBegin_s > MAX_GESTURE_DURATION_S))
    {
      if (gestureEnd_s - gestureBegin_s > MAX_GESTURE_DURATION_S)
      {
        // Take the initial time for the gesture begin.
        gestureBegin_s = openingBegin_s - 0.030;
      }
      hasError = true;
    }
  }
  
  // ****************************************************************
  // If the velic port is closed at the intended end of the opening,
  // move the end of the gesture rightwards in small steps, until the 
  // opening just happens.
  // ****************************************************************

  hasOpening = hasVelicOpening(gestureBegin_s, gestureEnd_s, openingEnd_s);
  hasError = false;

  while ((hasOpening == false) && (hasError == false))
  {
    // Move the end of the gesture 5 ms to the right.
    gestureEnd_s += TIME_STEP_S;
    hasOpening = hasVelicOpening(gestureBegin_s, gestureEnd_s, openingEnd_s);

    if (gestureEnd_s - gestureBegin_s > MAX_GESTURE_DURATION_S)
    {
      gestureEnd_s = gestureBegin_s + 0.010;
      hasError = true;
    }
  }

  // ****************************************************************
  // Put the new gesture into the tier.
  // ****************************************************************

  g.duration_s = gestureEnd_s - gestureBegin_s;
  g.dVal = TARGET_VELIC_OPENING;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = "";
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[VELIC_GESTURE].putGesture(g, gestureBegin_s);
}


// ****************************************************************************
/// Returns true, if the velic port is open at testTime_s, when a velic gesture
/// starting at gestureBegin_s and ending at gestureEnd_s is put into the
/// tier for velic gestures.
// ****************************************************************************

bool GesturalScore::hasVelicOpening(double gestureBegin_s, double gestureEnd_s, double testTime_s)
{
  const double VELIC_OPENING_THRESHOLD = 0.01;  // Minimally open port (1% of maximum value).

  static GestureSequence storedGestures;
  static Gesture g;
  static double tractParams[VocalTract::NUM_PARAMS];
  static double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  bool hasOpening = false;

  // Store the old state of the velic gesture tier.

  storedGestures = gestures[VELIC_GESTURE];

  // ****************************************************************
  // Put the gesture into the tier.
  // ****************************************************************

  g.duration_s = gestureEnd_s - gestureBegin_s;
  g.dVal = TARGET_VELIC_OPENING;
  g.neutral = false;
  g.slope = 0.0;
  g.sVal = "";
  g.tau_s = DEFAULT_TIME_CONSTANT_S;

  gestures[VELIC_GESTURE].putGesture(g, gestureBegin_s);

  // ****************************************************************
  // Calculate the parameter curves, the vocal tract parameters,
  // and test the parameter VO.
  // ****************************************************************

  calcCurves();
  getParams(testTime_s, tractParams, glottisParams);

  if (tractParams[VocalTract::VO] >= VELIC_OPENING_THRESHOLD)
  {
    hasOpening = true;
  }

  // Restore the old state of the gesture tier.
  gestures[VELIC_GESTURE] = storedGestures;

  return hasOpening;
}


// ****************************************************************************
/// Loads the gestural score from an XML-file.
// ****************************************************************************

bool GesturalScore::loadGesturesXml(const string &fileName, bool &allValuesInRange)
{
  // Assume here that all gesture values are in their valid ranges.
  allValuesInRange = true;

  // ****************************************************************
  // Load and parse the XML data.
  // ****************************************************************

  vector<XmlError> xmlErrors;
  XmlNode *rootNode = xmlParseFile(fileName, "gestural_score", &xmlErrors);
  if (rootNode == NULL)
  {
    xmlPrintErrors(xmlErrors);
    return false;
  }

  // ****************************************************************
  // Read the individual gesture sequences.
  // ****************************************************************

  XmlNode *node = NULL;
  int i, k;
  int N = rootNode->numChildElements("gesture_sequence");
  string st;
  bool theseValuesInRange = true;

  for (i=0; i < N; i++)
  {
    node = rootNode->getChildElement("gesture_sequence", i);
    st = node->getAttributeString("type");

    for (k=0; k < NUM_GESTURE_TYPES; k++)
    {
      if (gestures[k].abbr == st)
      {
        gestures[k].readFromXml(*node, theseValuesInRange);
        if (theseValuesInRange == false)
        {
          allValuesInRange = false;
          printf("Gesture values out of range for gesture type %d.\n", k);
        }
      }
    }
  }

  // Free the memory of the XML tree !
  delete rootNode;

  // Reset the synthesis after loading the score !
  resetSequence();

  return true;
}


// ****************************************************************************
/// Saves the gestural score in XML-format.
// ****************************************************************************

bool GesturalScore::saveGesturesXml(const string &fileName)
{
  ofstream os(fileName);
  int i;

  if (!os)
  {
    printf("Error: The file %s could not be opened!\n", fileName.c_str());
    return false;
  }

  // ****************************************************************
  // Open the <gestural-score> element and write the data.
  // ****************************************************************

  os << "<gestural_score>" << endl;
  
  for (i=0; i < NUM_GESTURE_TYPES; i++)
  {
    gestures[i].writeToXml(os, 2);
  }
  
  // ****************************************************************
  // Close the <gestural-score> element.
  // ****************************************************************

  os << "</gestural_score>" << endl;

  // Close the file

  os.close();

  return true;
}


// ****************************************************************************
/// Calculate the target sequences and the parameter curves from the gesture
/// sequences.
// ****************************************************************************

void GesturalScore::calcCurves()
{
  int i, k;
  double neutral;

  // ****************************************************************
  // Convert the gestures into target sequences.
  // ****************************************************************

  calcTractParamTargets();
  calcGlottisParamTargets();

  // ****************************************************************
  // Clear all curves with neutral values.
  // ****************************************************************

  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    neutral = vocalTract->param[i].neutral;
    for (k=0; k < MAX_CURVE_SAMPLES; k++)
    {
      tractParamCurve[i][k] = neutral;
    }
  }

  int numGlottisParams = (int)glottis->controlParam.size();
  for (i=0; i < numGlottisParams; i++)
  {
    neutral = glottis->controlParam[i].neutral;
    // For F0, initialize the neutral value in st.
    // They will be converted back into Hz below after filtering the st-targets.
    if (i == Glottis::FREQUENCY)
    {
      neutral = getF0_st(neutral);
    }

    for (k=0; k < MAX_CURVE_SAMPLES; k++)
    {
      glottisParamCurve[i][k] = neutral;
    }
  }

  // ****************************************************************
  // Convert the targets into time signals.
  // ****************************************************************

  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    calcParamCurve(tractParamTargets[i], tractParamCurve[i]);
  }

  numGlottisParams = (int)glottis->controlParam.size();
  for (i=0; i < numGlottisParams; i++)
  {
    calcParamCurve(glottisParamTargets[i], glottisParamCurve[i]);
  }

  // ****************************************************************
  // The targets and their slopes for F0 are expressed in semitones,
  // but we want the final parameter curve to be in Hz. So, transform
  // the values here.
  // ****************************************************************

  int numSamples = (int)glottisParamCurve[Glottis::FREQUENCY].size();
  for (i=0; i < numSamples; i++)
  {
    glottisParamCurve[Glottis::FREQUENCY][i] = getF0_Hz( glottisParamCurve[Glottis::FREQUENCY][i] );
  }

}


// ****************************************************************************
/// Return the parameter values of the parameter curves for the glottis and the
/// vocal tract.
/// Thereby, interpolate between adjacent curve sample points.
// ****************************************************************************

void GesturalScore::getParams(double pos_s, double *vocalTractParams, double *glottisParams)
{
  const double SAMPLING_PERIOD_S = 1.0 / (double)CURVE_SAMPLING_RATE;
  int index = (int)(pos_s * (double)CURVE_SAMPLING_RATE);
  int numSamples = (int)tractParamCurve[0].size();   // numSamples should be equal for all parameter curves !!!
  double s = (pos_s - index*SAMPLING_PERIOD_S) / SAMPLING_PERIOD_S;
  double s1 = 1.0 - s;
  int i;

  // For safety !
  if (index < 0)
  {
    index = 0;
  }
  
  if (index > MAX_CURVE_SAMPLES-2)
  {
    index = MAX_CURVE_SAMPLES-2;
    s = 0.0;
    s1 = 1.0;
  }

  if (index >= numSamples-1)
  {
    index = numSamples - 2;
    s = 1.0;
    s1 = 0.0;
  }

  if (vocalTractParams != NULL)
  {
    for (i=0; i < VocalTract::NUM_PARAMS; i++)
    {
      vocalTractParams[i] = s1*tractParamCurve[i][index] + s*tractParamCurve[i][index + 1];
    }
  }

  if (glottisParams != NULL)
  {
    int numGlottisParams = (int)glottis->controlParam.size();
    for (i=0; i < numGlottisParams; i++)
    {
      glottisParams[i] = s1*glottisParamCurve[i][index] + s*glottisParamCurve[i][index + 1];
    }
  }

}


// ****************************************************************************
/// Returns the duration of the gestural score in seconds. This is the duration
/// of the longest gesture sequence in the score.
// ****************************************************************************

double GesturalScore::getScoreDuration_s()
{
  double duration_s = 0.0;
  double sequenceDur_s = 0.0;
  int i;

  for (i=0; i < NUM_GESTURE_TYPES; i++)
  {
    sequenceDur_s = gestures[i].getDuration_s();
    if (sequenceDur_s > duration_s)
    {
      duration_s = sequenceDur_s;
    }
  }

  return duration_s;
}


// ****************************************************************************
/// Moves the whole F0 curve up or down by the given offset in semitones.
// ****************************************************************************

void GesturalScore::changeF0Offset(double deltaF0_st)
{
  printf("Changing F0 offset by %2.2f st.\n", deltaF0_st);

  GestureSequence *sequence = &gestures[F0_GESTURE];
  Gesture *g;
  int N = sequence->numGestures();
  int i;

  for (i = 0; i < N; i++)
  {
    g = sequence->getGesture(i);
    g->dVal += deltaF0_st;

    if (g->dVal < sequence->minValue)
    {
      g->dVal = sequence->minValue;
      printf("F0 target value has been limited to %2.4f.\n", sequence->minValue);
    }

    if (g->dVal > sequence->maxValue)
    {
      g->dVal = sequence->maxValue;
      printf("F0 target value has been limited to %2.4f.\n", sequence->maxValue);
    }
  }

  // Important:
  calcCurves();
}


// ****************************************************************************
/// Changes the range of the F0 targets by the given factor.
/// Afterwards the whole contour is moved up ur down such that the average F0
/// in semitones is the same as before.
// ****************************************************************************

void GesturalScore::changeF0Range(double factor)
{
  // Get the current (old) F0 statistics.
  double oldF0Mean_st;
  double oldF0Sd_st;
  double oldF0Mean_Hz;
  double oldF0Sd_Hz;

  getF0Statistic(oldF0Mean_st, oldF0Sd_st, oldF0Mean_Hz, oldF0Sd_Hz);

  // ****************************************************************
  // Scale all F0 target values and slopes up by the given factor.
  // ****************************************************************

  printf("Scaling F0 range by the factor %2.2f.\n", factor);

  GestureSequence *sequence = &gestures[F0_GESTURE];
  Gesture *g;
  int N = sequence->numGestures();
  int i;

  for (i = 0; i < N; i++)
  {
    g = sequence->getGesture(i);
    g->dVal *= factor;
    g->slope *= factor;

    // Don't check the boundaries of the F0 values, because they
    // may temporarily exceed the boundaries, because the whole
    // F0 contour is shifted (corrected) further below!

    if (g->slope < sequence->minSlope)
    {
      g->slope = sequence->minSlope;
      printf("F0 slope has been limited to %2.4f.\n", sequence->minSlope);
    }

    if (g->slope > sequence->maxSlope)
    {
      g->slope = sequence->maxSlope;
      printf("F0 slope has been limited to %2.4f.\n", sequence->maxSlope);
    }
  }

  // Get the new F0 statistics.
  double newF0Mean_st;
  double newF0Sd_st;
  double newF0Mean_Hz;
  double newF0Sd_Hz;

  getF0Statistic(newF0Mean_st, newF0Sd_st, newF0Mean_Hz, newF0Sd_Hz);

  // ****************************************************************
  // Correct the induced change in the mean F0.
  // ****************************************************************

  double f0MeanDiff_st = newF0Mean_st - oldF0Mean_st;
  printf("The mean F0 was changed by %2.2f st and will hence be corrected by %2.2f st.\n",
    f0MeanDiff_st,
    -f0MeanDiff_st);

  changeF0Offset(-f0MeanDiff_st);

  // Important:
  calcCurves();
}


// ****************************************************************************
/// Changes the slope of all F0 gestures by the given delta value.
// ****************************************************************************

void GesturalScore::changeF0TargetSlope(double deltaSlope_st_s)
{
  // Get the current (old) F0 statistics.
  double oldF0Mean_st;
  double oldF0Sd_st;
  double oldF0Mean_Hz;
  double oldF0Sd_Hz;

  getF0Statistic(oldF0Mean_st, oldF0Sd_st, oldF0Mean_Hz, oldF0Sd_Hz);

  // ****************************************************************
  // Changing all F0 target slopes by adding the given number of semitones.
  // ****************************************************************

  printf("Changing all F0 target slopes by adding %2.2f semitones.\n", deltaSlope_st_s);

  GestureSequence *sequence = &gestures[F0_GESTURE];
  Gesture *g;
  int N = sequence->numGestures();
  int i;

  for (i = 0; i < N; i++)
  {
    g = sequence->getGesture(i);
    g->slope += deltaSlope_st_s;

    // Don't check the boundaries of the F0 values, because they
    // may temporarily exceed the boundaries, because the whole
    // F0 contour is shifted (corrected) further below!

    if (g->slope < sequence->minSlope)
    {
      g->slope = sequence->minSlope;
      printf("F0 slope has been limited to %2.4f.\n", sequence->minSlope);
    }

    if (g->slope > sequence->maxSlope)
    {
      g->slope = sequence->maxSlope;
      printf("F0 slope has been limited to %2.4f.\n", sequence->maxSlope);
    }
  }

  // Get the new F0 statistics.
  double newF0Mean_st;
  double newF0Sd_st;
  double newF0Mean_Hz;
  double newF0Sd_Hz;

  getF0Statistic(newF0Mean_st, newF0Sd_st, newF0Mean_Hz, newF0Sd_Hz);

  // ****************************************************************
  // Correct the induced change in the mean F0.
  // ****************************************************************

  double f0MeanDiff_st = newF0Mean_st - oldF0Mean_st;
  printf("The mean F0 was changed by %2.2f st and will hence be corrected by %2.2f st.\n",
    f0MeanDiff_st,
    -f0MeanDiff_st);

  changeF0Offset(-f0MeanDiff_st);

  // Important:
  calcCurves();
}


// ****************************************************************************
/// Substitutes all glottal shape gestures with the name oldShapeName by
/// newShapeName, e.g. all breathy phonation gestures by model phonation gestures.
// ****************************************************************************

void GesturalScore::substituteGlottalShapes(const string &oldShapeName, const string &newShapeName)
{
  printf("Substituting glottal shape '%s' by '%s'.\n", oldShapeName.c_str(), newShapeName.c_str());

  GestureSequence *sequence = &gestures[GLOTTAL_SHAPE_GESTURE];
  Gesture *g;
  int N = sequence->numGestures();
  int i;

  for (i = 0; i < N; i++)
  {
    g = sequence->getGesture(i);
    if (g->sVal == oldShapeName)
    {
      g->sVal = newShapeName;
    }
  }

  // Important:
  calcCurves();
}


// ****************************************************************************
/// Scale the subglottal pressure by the given factor.
// ****************************************************************************

void GesturalScore::changeSubglottalPressure(double factor)
{
  printf("Changing subglottal pressure by the factor %2.2f.\n", factor);

  GestureSequence *sequence = &gestures[PRESSURE_GESTURE];
  Gesture *g;
  int N = sequence->numGestures();
  int i;

  for (i = 0; i < N; i++)
  {
    g = sequence->getGesture(i);
    g->dVal *= factor;

    if (g->dVal < sequence->minValue)
    {
      g->dVal = sequence->minValue;
      printf("Pressure has been limited to %d dPa.\n", (int)sequence->minValue);
    }

    if (g->dVal > sequence->maxValue)
    {
      g->dVal = sequence->maxValue;
      printf("Pressure has been limited to %d dPa.\n", (int)sequence->maxValue);
    }
  }

  // Important:
  calcCurves();
}


// ****************************************************************************
/// Scale the length of all gestures with the given factor, such that the 
/// utterance is spoken slower or faster. The slopes of F0 gestures are also
/// changed accordingly.
// ****************************************************************************

void GesturalScore::changeDuration(double factor)
{
  if (factor > 4.0)
  {
    factor = 4.0;
    printf("Factor for change of gestural score duration has been limited to 4.0.\n");
  }

  if (factor < 0.25)
  {
    factor = 0.25;
    printf("Factor for change of gestural score duration has been limited to 0.25.\n");
  }

  printf("Changing gestural score duration by the factor %2.2f.\n", factor);

  int i, k;
  GestureSequence *sequence = NULL;
  Gesture *gesture = NULL;
  int numGestures = 0;

  for (i = 0; i < NUM_GESTURE_TYPES; i++)
  {
    sequence = &gestures[i];
    numGestures = sequence->numGestures();

    for (k = 0; k < numGestures; k++)
    {
      gesture = sequence->getGesture(k);
      gesture->duration_s *= factor;

      if (i == F0_GESTURE)
      {
        gesture->slope /= factor;

        if (gesture->slope < sequence->minSlope)
        {
          gesture->slope = sequence->minSlope;
          printf("F0 target slope has been limited to %2.4f.\n", sequence->minSlope);
        }

        if (gesture->slope > sequence->maxSlope)
        {
          gesture->slope = sequence->maxSlope;
          printf("F0 target slope has been limited to %2.4f.\n", sequence->maxSlope);
        }
      }
    }
  }

  calcCurves();
}


// ****************************************************************************
/// Change the time constant of all gestures by the given factor.
// ****************************************************************************

void GesturalScore::changeTimeConstants(double factor)
{
  printf("Changing all time constants in the gestural score by the factor %2.2f.\n", factor);

  int i, k;
  GestureSequence *sequence = NULL;
  Gesture *g = NULL;
  int numGestures = 0;

  for (i = 0; i < NUM_GESTURE_TYPES; i++)
  {
    sequence = &gestures[i];
    numGestures = sequence->numGestures();

    for (k = 0; k < numGestures; k++)
    {
      g = sequence->getGesture(k);
      g->tau_s *= factor;

      if (g->tau_s < sequence->minTau_s)
      {
        g->tau_s = sequence->minTau_s;
        printf("Time constant has been limited to %2.4f.\n", sequence->minTau_s);
      }

      if (g->tau_s > sequence->maxTau_s)
      {
        g->tau_s = sequence->maxTau_s;
        printf("Time constant has been limited to %2.4f.\n", sequence->maxTau_s);
      }
    }
  }

  calcCurves();
}


// ****************************************************************************
/// Obtains the given reference parameters form the model F0 curve.
// ****************************************************************************

void GesturalScore::getF0Statistic(double &f0Mean_st, double &f0Sd_st, double &f0Mean_Hz, double &f0Sd_Hz)
{
  int i;
  double d;

  calcCurves();

  // ****************************************************************
  // After calcCurves(), the values in 
  // glottisParamCurve[Glottis::FREQUENCY][i] are in Hz !
  // ****************************************************************

  int numSamples = (int)(getScoreDuration_s() * CURVE_SAMPLING_RATE) - 1;

  if (numSamples < 1)
  {
    f0Mean_Hz = 0.0;
    f0Sd_Hz = 0.0;
    f0Mean_st = 0.0;
    f0Sd_st = 0.0;
    return;
  }
//  printf("Calculating gestural score F0 statistics for a duration of %2.4f s.\n", (double)numSamples / (double)CURVE_SAMPLING_RATE);

  f0Mean_Hz = 0.0;
  f0Sd_Hz = 0.0;
  for (i = 0; i < numSamples; i++)
  {
    f0Mean_Hz += glottisParamCurve[Glottis::FREQUENCY][i];
  }

  f0Mean_Hz /= (double)numSamples;

  for (i = 0; i < numSamples; i++)
  {
    d = glottisParamCurve[Glottis::FREQUENCY][i] - f0Mean_Hz;
    f0Sd_Hz += d*d;
  }
  f0Sd_Hz /= (double)numSamples;
  f0Sd_Hz = sqrt(f0Sd_Hz);  

  // ****************************************************************
  // Transform the F0 values to semitones and then calculate mean
  // and standard deviation.
  // ****************************************************************

  for (i = 0; i < numSamples; i++)
  {
    glottisParamCurve[Glottis::FREQUENCY][i] = getF0_st(glottisParamCurve[Glottis::FREQUENCY][i]);
  }

  f0Mean_st = 0.0;
  f0Sd_st = 0.0;
  for (i = 0; i < numSamples; i++)
  {
    f0Mean_st += glottisParamCurve[Glottis::FREQUENCY][i];
  }

  f0Mean_st /= (double)numSamples;

  for (i = 0; i < numSamples; i++)
  {
    d = glottisParamCurve[Glottis::FREQUENCY][i] - f0Mean_st;
    f0Sd_st += d*d;
  }
  f0Sd_st /= (double)numSamples;
  f0Sd_st = sqrt(f0Sd_st);

  // Output the results.

//  printf("F0 mean (SD) with a Hz scale: %2.2f (%2.2f) Hz\n", f0Mean_Hz, f0Sd_Hz);
//  printf("F0 mean (SD) with a semitone scale: %2.2f (%2.2f) st\n", f0Mean_st, f0Sd_st);

  // Important: Calculate the F0 samples in Hz again!
  calcCurves();
}


// ****************************************************************************
/// Calculates the pseudo-inverse matrix of the matrix M with numRows rows and
/// two columns. The result is a two-dimensional array that is pointed to by
/// invM.
// ****************************************************************************

void GesturalScore::getPseudoInverseNx2(double M[][2], int numRows, double *invM)
{
  int x, y, k;

  // ****************************************************************
  // Calculate M^T * M (which is a 2x2 matrix).
  // ****************************************************************

  double A[2][2];
  double sum = 0.0;

  for (x = 0; x < 2; x++)
  {
    for (y = 0; y < 2; y++)
    {
      sum = 0.0;
      for (k = 0; k < numRows; k++)
      {
        sum += M[k][y] * M[k][x];
      }
      A[x][y] = sum;
    }
  }

  // ****************************************************************
  // Invert A to obtain B.
  // ****************************************************************

  double B[2][2];
  const double EPSILON = 0.000000001;

  double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  if (fabs(det) < EPSILON)
  {
    det = EPSILON;
  }
  B[0][0] = A[1][1] / det;
  B[0][1] = -A[0][1] / det;
  B[1][0] = -A[1][0] / det;
  B[1][1] = A[0][0] / det;

  // ****************************************************************
  // Calculate B * M^T to obtain the pseudo-inverse of M.
  // ****************************************************************

  for (x = 0; x < numRows; x++)
  {
    for (y = 0; y < 2; y++)
    {
      invM[y*numRows + x] = B[y][0] * M[x][0] + B[y][1] * M[x][1];
    }
  }

}



// ****************************************************************************
/// Maps the given vocal tract shape into the subspace spanned by the three
/// corner vowels /a/, /i/, and /u/. The tongue shape and the lip shape are 
/// mapped individually so that lip rounding is independent from tongue position.
/// The coordinates in vowel subspace are alpha and beta.
/// x = x_a + alpha*(x_i - x_a) + beta*(x_u - x_a).
/// The mapping is done using the pseudeinverse of the coefficient matrix.
// ****************************************************************************

void GesturalScore::mapToVowelSubspace(VocalTract *vt, double *vocalTractParams, 
  double &alphaTongue, double &betaTongue, double &alphaLips, double &betaLips)
{
  int i;

  alphaTongue = 0.0;
  betaTongue = 0.0;
  alphaLips = 0.0;
  betaLips = 0.0;

  // Get the indices of the vocal tract shapes for the vowels /i/, /a/, and /u/.

  int shapeIndexA = vt->getShapeIndex("a");
  int shapeIndexI = vt->getShapeIndex("i");
  int shapeIndexU = vt->getShapeIndex("u");

  // Return if one of the shapes was not found.

  if ((shapeIndexA == -1) || (shapeIndexI == -1) || (shapeIndexU == -1))
  {
    return;
  }

  // ****************************************************************
  // Map the TONGUE SHAPE of the given shape to the subspace
  // defined by the shapes for /a/, /i/, and /u/.
  // ****************************************************************

  // The list of vocal tract parameters to map the tongue shape.

  const int NUM_TONGUE_PARAMS = 17;
  const int TONGUE_PARAM_INDEX[NUM_TONGUE_PARAMS] =
  {
    VocalTract::HX, 
    VocalTract::HY, 
    VocalTract::JX, 
    VocalTract::JA, 
    VocalTract::VS, 
    VocalTract::VO, 
    VocalTract::TCX, 
    VocalTract::TCY, 
    VocalTract::TTX, 
    VocalTract::TTY,
    VocalTract::TBX, 
    VocalTract::TBY, 
    VocalTract::TRX, 
    VocalTract::TRY,
    VocalTract::TS1, 
    VocalTract::TS2, 
    VocalTract::TS3
  };

  // Calculate matrix and solution vector for the first LSE.
  double M[NUM_TONGUE_PARAMS][2];     // Mapping matrix
  double v[NUM_TONGUE_PARAMS];        // Solution vector

  for (i=0; i < NUM_TONGUE_PARAMS; i++)
  {
    // Matrix column for the vector x_i - x_a.
    M[i][0] = 
      vt->shapes[shapeIndexI].param[ TONGUE_PARAM_INDEX[i] ] - 
      vt->shapes[shapeIndexA].param[ TONGUE_PARAM_INDEX[i] ];

    // Matrix column for the vector x_u - x_a.
    M[i][1] = 
      vt->shapes[shapeIndexU].param[ TONGUE_PARAM_INDEX[i] ] - 
      vt->shapes[shapeIndexA].param[ TONGUE_PARAM_INDEX[i] ];

    // Solution vector = x_given - x_a
    v[i] = 
      vocalTractParams[ TONGUE_PARAM_INDEX[i] ] - 
      vt->shapes[shapeIndexA].param[ TONGUE_PARAM_INDEX[i] ];
  }

  // The overdetermined system of equations is now M*(alpha beta)^T = b.
  // To get the best estimate of alpha, beta in the least square sense,
  // calculate (alpha beta)^T = inv(M)*b.

  double invM[2][NUM_TONGUE_PARAMS];
  getPseudoInverseNx2(M, NUM_TONGUE_PARAMS, &invM[0][0]);

  alphaTongue = 0.0;
  betaTongue = 0.0;
  for (i = 0; i < NUM_TONGUE_PARAMS; i++)
  {
    alphaTongue += invM[0][i] * v[i];
    betaTongue += invM[1][i] * v[i];
  }

  
  // ****************************************************************
  // Map the LIP SHAPE of the given shape to the subspace
  // defined by the shapes for /a/, /i/, and /u/.
  // ****************************************************************

  // Let the coefficient matrix be (a b; c d) and the solution vector
  // is (x y).

  double a = vt->shapes[shapeIndexI].param[ VocalTract::LP ] - 
             vt->shapes[shapeIndexA].param[ VocalTract::LP ];

  double b = vt->shapes[shapeIndexU].param[ VocalTract::LP ] - 
             vt->shapes[shapeIndexA].param[ VocalTract::LP ];

  double c = vt->shapes[shapeIndexI].param[ VocalTract::LD ] - 
             vt->shapes[shapeIndexA].param[ VocalTract::LD ];

  double d = vt->shapes[shapeIndexU].param[ VocalTract::LD ] - 
             vt->shapes[shapeIndexA].param[ VocalTract::LD ];

  double x = vocalTractParams[ VocalTract::LP ] - vt->shapes[shapeIndexA].param[ VocalTract::LP ];
  double y = vocalTractParams[ VocalTract::LD ] - vt->shapes[shapeIndexA].param[ VocalTract::LD ];

  const double EPSILON = 0.000000001;
  double det = a*d - b*c;
  if (fabs(det) < EPSILON)
  {
    det = EPSILON;
  }

  alphaLips = (x*d - b*y) / det;
  betaLips = (a*y - x*c) / det;
}


// ****************************************************************************
/// Limits the coordinates (alpha, beta) of a vowel on the a-i-u-subspace into 
/// the region that is defined by the lines between /a/, /i/, and /u/.
/// That means, we require that 0<=alpha<=1; 0<=beta<=1; alpha+beta<=1.
// ****************************************************************************

void GesturalScore::limitVowelSubspaceCoord(double &alphaTongue, double &betaTongue, double &alphaLips, double &betaLips)
{
  // ****************************************************************
  // For the tongue.
  // ****************************************************************

  if (alphaTongue < 0.0)
  {
    alphaTongue = 0.0;
  }

  if (alphaTongue > 1.0)
  {
    alphaTongue = 1.0;
  }

  if (betaTongue < 0.0)
  {
    betaTongue = 0.0;
  }

  if (betaTongue > 1.0)
  {
    betaTongue = 1.0;
  }

  if (alphaTongue + betaTongue > 1.0)
  {
    double d = alphaTongue + betaTongue - 1.0;
    alphaTongue-= 0.5*d;
    betaTongue-= 0.5*d;
  }

  // ****************************************************************
  // For the lips.
  // ****************************************************************

  if (alphaLips < 0.0)
  {
    alphaLips = 0.0;
  }

  if (alphaLips > 1.0)
  {
    alphaLips = 1.0;
  }

  if (betaLips < 0.0)
  {
    betaLips = 0.0;
  }

  if (betaLips > 1.0)
  {
    betaLips = 1.0;
  }

  if (alphaLips + betaLips > 1.0)
  {
    double d = alphaLips + betaLips - 1.0;
    alphaLips-= 0.5*d;
    betaLips-= 0.5*d;
  }
}


// ****************************************************************************
/// Interpolates the consonantal vocal tract shape with the given name (e.g., 
/// "tt-alveolar-stop") between the three prototypes for this consonant in the
/// context of the vowels /a/, /i/, and /u/. The tongue and lips are 
/// interpolated independently with the corresponding alphas and betas as
/// x = x_a + alpha*(x_i - x_a) + beta*(x_u - x_a).
// ****************************************************************************

bool GesturalScore::getContextDependentConsonant(VocalTract *vt, const char *name, 
  double alphaTongue, double betaTongue, double alphaLips, double betaLips, 
  double *consonantParams)
{
  int i;

  std::string n(name);
  int consonantIndexA = vt->getShapeIndex( n + "(a)" );
  int consonantIndexI = vt->getShapeIndex( n + "(i)" );
  int consonantIndexU = vt->getShapeIndex( n + "(u)" );
  
  if ((consonantIndexA == -1) || (consonantIndexI == -1) || (consonantIndexU == -1))
  {
    // Set the neutral vocal tract.
    for (i=0; i < VocalTract::NUM_PARAMS; i++)
    {
      consonantParams[i] = vt->param[i].neutral;
    }
    return false;
  }

  // ****************************************************************

  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    consonantParams[i] = vt->shapes[consonantIndexA].param[i] +
      alphaTongue*(vt->shapes[consonantIndexI].param[i] - vt->shapes[consonantIndexA].param[i]) +
      betaTongue*(vt->shapes[consonantIndexU].param[i] - vt->shapes[consonantIndexA].param[i]);
  }

  // The two lip parameters are interpolated independently with 
  // alphaLips and betaLips.

  i = VocalTract::LP;
  consonantParams[i] = vt->shapes[consonantIndexA].param[i] +
    alphaLips*(vt->shapes[consonantIndexI].param[i] - vt->shapes[consonantIndexA].param[i]) +
    betaLips*(vt->shapes[consonantIndexU].param[i] - vt->shapes[consonantIndexA].param[i]);

  i = VocalTract::LD;
  consonantParams[i] = vt->shapes[consonantIndexA].param[i] +
    alphaLips*(vt->shapes[consonantIndexI].param[i] - vt->shapes[consonantIndexA].param[i]) +
    betaLips*(vt->shapes[consonantIndexU].param[i] - vt->shapes[consonantIndexA].param[i]);

  return true;
}


// ****************************************************************************
/// Returns the vocal tract parameters of the consonant consonantName 
/// (e.g., tt-alveolar-stop) (co)articulated in the context of the vowel 
/// vowelName.
// ****************************************************************************

bool GesturalScore::getContextDependentConsonant(VocalTract *vt,
  const char *consonantName, const char *contextVowelName, double *consonantParams)
{
  int k;

  // ****************************************************************
  // Obtain the vocal tract parameters of the context vowel.
  // ****************************************************************

  double contextVowelParams[VocalTract::NUM_PARAMS];

  int shapeIndex = vt->getShapeIndex(string(contextVowelName));

  // Try to get the index of Schwa.
  if (shapeIndex == -1)
  {
    shapeIndex = vt->getShapeIndex("@");
  }

  if (shapeIndex == -1)
  {
    for (k = 0; k < VocalTract::NUM_PARAMS; k++)
    {
      contextVowelParams[k] = vt->param[k].neutral;
    }
  }
  else
  {
    for (k = 0; k < VocalTract::NUM_PARAMS; k++)
    {
      contextVowelParams[k] = vt->shapes[shapeIndex].param[k];
    }
  }

  // ****************************************************************
  // Get the coordinates of the context vowel in the 2D-vowel subspace.
  // ****************************************************************

  double alphaTongue;
  double betaTongue;
  double alphaLips;
  double betaLips;

  mapToVowelSubspace(vt, contextVowelParams, alphaTongue, betaTongue, alphaLips, betaLips);
  limitVowelSubspaceCoord(alphaTongue, betaTongue, alphaLips, betaLips);

  // ****************************************************************
  // Get the context-dependent consonant parameters.
  // ****************************************************************

  bool ok = getContextDependentConsonant(vt, consonantName,
    alphaTongue, betaTongue, alphaLips, betaLips,
    consonantParams);

  return ok;
}


// ****************************************************************************
/// Returns the cross-sectional area of the constriction in the vocal tract
/// for the given vocal tract parameters at the region of the area function
/// that is typical for the given (consonantal) gesture type.
/// Supported gesture types are: 
/// LIP_GESTURE, TONGUE_TIP_GESTURE, TONGUE_BODY_GESTURE.
// ****************************************************************************

double GesturalScore::getConstrictionArea_cm2(VocalTract *vt,
  double *vocalTractParams, GestureType gestureType)
{
  int i;

  // ****************************************************************
  // Calculate the vocal tract shape for the given parameters.
  // ****************************************************************

  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vt->param[i].x = vocalTractParams[i];
  }
  vt->calculateAll();

  // ****************************************************************
  // Determine the minimal cross-sectional area depending on the 
  // gesture type.
  // ****************************************************************

  double constrictionArea_cm2 = 1000000.0;   // Very big value by default.
  VocalTract::TubeSection *ts = NULL;

  if (gestureType == LIP_GESTURE)
  {
    for (i = 0; i < VocalTract::NUM_TUBE_SECTIONS; i++)
    {
      ts = &vt->tubeSection[i];
      if ((ts->articulator == Tube::LOWER_LIP) && (ts->area < constrictionArea_cm2))
      {
        constrictionArea_cm2 = ts->area;
      }
    }
  }
  else

  // We look for a constriction formed with the tongue body or tongue tip.
  {
    // Determine the first section of the "tongue tip part".
    // Therefore, find the last tongue section first.

    int lastTongueSection = 0;
    for (i = 1; i < VocalTract::NUM_TUBE_SECTIONS; i++)
    {
      ts = &vt->tubeSection[i];
      if (ts->articulator == Tube::TONGUE)
      {
        lastTongueSection = i;
      }
    }

    int firstTongueTipSection = lastTongueSection;
    double summedLength_cm = vt->tubeSection[firstTongueTipSection].length;

    // Count so many previous sections to the tongue tip that the
    // tongue tip is at least 2 cm long.

    while ((firstTongueTipSection > 0) && (summedLength_cm < 2.0))
    {
      firstTongueTipSection--;
      summedLength_cm += vt->tubeSection[firstTongueTipSection].length;
    }

    // Determine the minimum area in the respective regions.

    if (gestureType == TONGUE_TIP_GESTURE)
    {
      for (i = firstTongueTipSection; i < VocalTract::NUM_TUBE_SECTIONS; i++)
      {
        ts = &vt->tubeSection[i];
        if ((ts->articulator == Tube::TONGUE) && (ts->area < constrictionArea_cm2))
        {
          constrictionArea_cm2 = ts->area;
        }
      }
    }
    else
    if (gestureType == TONGUE_BODY_GESTURE)
    {
      for (i = 0; i < firstTongueTipSection; i++)
      {
        ts = &vt->tubeSection[i];
        if ((ts->articulator == Tube::TONGUE) && (ts->area < constrictionArea_cm2))
        {
          constrictionArea_cm2 = ts->area;
        }
      }
    }

  }

  return constrictionArea_cm2;
}


// ****************************************************************************
// Returns the F0 given by freq_Hz in semitones relative to REFERENCE_FREQUENCY.
// ****************************************************************************

double GesturalScore::getF0_st(double freq_Hz)
{
  if (freq_Hz < 1.0) { freq_Hz = 1.0; }
  return 12.0*log(freq_Hz / REFERENCE_FREQUENCY) / log(2.0);
}


// ****************************************************************************
// Transforms a frequency from semitones in Hz.
// ****************************************************************************

double GesturalScore::getF0_Hz(double freq_st)
{
  return REFERENCE_FREQUENCY*pow(2, freq_st / 12.0);
}


// ****************************************************************************
/// Returns the complete tube geometry for the current position.
// ****************************************************************************

void GesturalScore::getTube(Tube &tube)
{
  int i;
  double pos_s = (double)pos / (double)SAMPLING_RATE;

  // ****************************************************************
  // Get a new tube shape from the vocal tract model every 2.5 ms.
  // ****************************************************************

  const double CURVE_SAMPLING_PERIOD = 1.0 / (double)CURVE_SAMPLING_RATE;
  int neededLeftIndex;
  double ratio, ratio1;

  neededLeftIndex = (int)(pos_s * (double)CURVE_SAMPLING_RATE);
  if (neededLeftIndex > MAX_CURVE_SAMPLES-2)
  {
    neededLeftIndex = MAX_CURVE_SAMPLES-2;
  }

  ratio = (pos_s - neededLeftIndex*CURVE_SAMPLING_PERIOD) / CURVE_SAMPLING_PERIOD;
  ratio1 = 1.0 - ratio;

  // ****************************************************************
  // The actual tube is interpolated between the "left tube" and 
  // the "right tube". Make sure that both these tubes have the
  // correct shape.
  // ****************************************************************

  if (leftTubeIndex != neededLeftIndex)
  {
    if ((leftTubeIndex != -1) && (neededLeftIndex == leftTubeIndex + 1))
    {
      // Exchange pointers to the tubes -> make the old right tube the new left tube.
      Tube *dummy = leftTube;
      leftTube = rightTube;
      rightTube = dummy;

      // Get the new shape for the right tube.
      for (i=0; i < VocalTract::NUM_PARAMS; i++)
      {
        vocalTract->param[i].x = tractParamCurve[i][neededLeftIndex + 1];
      }
      vocalTract->calculateAll();
      vocalTract->getTube(rightTube);

      leftTubeIndex = neededLeftIndex;
    }
    else
    {
      // Get the left tube
      for (i=0; i < VocalTract::NUM_PARAMS; i++)
      {
        vocalTract->param[i].x = tractParamCurve[i][neededLeftIndex];
      }
      vocalTract->calculateAll();
      vocalTract->getTube(leftTube);

      // Get the right tube
      for (i=0; i < VocalTract::NUM_PARAMS; i++)
      {
        vocalTract->param[i].x = tractParamCurve[i][neededLeftIndex + 1];
      }
      vocalTract->calculateAll();
      vocalTract->getTube(rightTube);

      leftTubeIndex = neededLeftIndex;
    }
  }

  // ****************************************************************
  // Interpolate the tube.
  // ****************************************************************

  tube.interpolate(leftTube, rightTube, ratio);

  // ****************************************************************
  // Get the current glottis params and calculate the glottis 
  // geometry.
  // ****************************************************************

  double length_cm[Tube::NUM_GLOTTIS_SECTIONS];
  double area_cm2[Tube::NUM_GLOTTIS_SECTIONS];
  int numGlottisParams = (int)glottis->controlParam.size();
  double glottisParams[128];

  getParams(pos_s, NULL, glottisParams);

  for (i=0; i < numGlottisParams; i++)
  {
    glottis->controlParam[i].x = glottisParams[i];
  }
  glottis->calcGeometry();
  glottis->getTubeData(length_cm, area_cm2);

  tube.setGlottisGeometry(length_cm, area_cm2);
  tube.setAspirationStrength( glottis->getAspirationStrength_dB() );
}


// ****************************************************************************
/// This synthesis mode has only a pressure source.
// ****************************************************************************

void GesturalScore::getFlowSource(double &flow_cm3_s, int &section)
{
  flow_cm3_s = 0.0;
  section = -1;
}


// ****************************************************************************
// ****************************************************************************

void GesturalScore::getPressureSource(double &pressure_dPa, int &section)
{
  section = Tube::FIRST_TRACHEA_SECTION;
  pressure_dPa = glottis->controlParam[Glottis::PRESSURE].x;
}


// ****************************************************************************
/// Reset the state of the tube sequence.
// ****************************************************************************

void GesturalScore::resetSequence()
{
  pos = 0;
  leftTubeIndex = -1;
  calcCurves();
  glottis->resetMotion();
}


// ****************************************************************************
// ****************************************************************************

void GesturalScore::incPos(const double pressure_dPa[])
{
  // Increment the time/sample number
  glottis->incTime(1.0/(double)SAMPLING_RATE, pressure_dPa);
  pos++;
}


// ****************************************************************************
/// Returns the duration of the gestural score in sample points, i. e. the 
/// duration of the longest gesture sequence.
// ****************************************************************************

int GesturalScore::getDuration_pt()
{
  int i;
  double maxDur_s = 0.0;
  double dur_s;

  for (i=0; i < NUM_GESTURE_TYPES; i++)
  {
    dur_s = gestures[i].getDuration_s();
    if (dur_s > maxDur_s)
    {
      maxDur_s = dur_s;
    }
  }

  return (int)(maxDur_s*SAMPLING_RATE);
}


// ****************************************************************************
// ****************************************************************************

int GesturalScore::getPos_pt()
{
  return pos;
}

// ****************************************************************************
// ****************************************************************************

void GesturalScore::calcTractParamTargets()
{
  const int MAX_SLICES = 1024;
  const int NUM_TRACT_GESTURE_TYPES = 5;

  Gesture slice[MAX_SLICES][NUM_TRACT_GESTURE_TYPES];
  int i, k, m;
  double pos_s = 0.0;
  double nearestPos_s;
  double sliceDuration_s;
  int index[NUM_TRACT_GESTURE_TYPES] = { 0, 0, 0, 0, 0 };
  double gesturePos_s[NUM_TRACT_GESTURE_TYPES] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  GestureSequence *sequence[NUM_TRACT_GESTURE_TYPES];
  int numSlices = 0;
  Gesture g;

  for (i=0; i < NUM_TRACT_GESTURE_TYPES; i++)
  {
    sequence[i] = &gestures[i];
  }

  // ****************************************************************
  // Split the gestures in smaller pieces. Each boundary between two 
  // gestures in the original score is introduced in all other
  // gesture rows.
  // ****************************************************************

  while ((index[0] < sequence[0]->numGestures()) ||
    (index[1] < sequence[1]->numGestures()) ||
    (index[2] < sequence[2]->numGestures()) ||
    (index[3] < sequence[3]->numGestures()) ||
    (index[4] < sequence[4]->numGestures()))
  {
    // **************************************************************
    // Find the next nearest border from position pos_s.
    // **************************************************************

    k = -1;
    nearestPos_s = 1000000.0;
    
    for (i=0; i < NUM_TRACT_GESTURE_TYPES; i++)
    {
      if (index[i] < sequence[i]->numGestures())
      {
        if (gesturePos_s[i] + sequence[i]->getGesture(index[i])->duration_s < nearestPos_s)
        {
          nearestPos_s = gesturePos_s[i] + sequence[i]->getGesture(index[i])->duration_s;
          k = i;
        }
      }
    }

    // **************************************************************
    // Add a gesture slices from pos_s to the end of the found 
    // gesture.
    // **************************************************************

    if (k != -1)
    {
      if (nearestPos_s > pos_s)
      {
        sliceDuration_s = nearestPos_s - pos_s;
        for (i=0; i < NUM_TRACT_GESTURE_TYPES; i++)
        {
          if (index[i] < sequence[i]->numGestures())
          {
            g = *sequence[i]->getGesture(index[i]);
          }
          else
          {
            // Create a neutral gesture
            g.dVal = sequence[i]->minValue;
            g.tau_s = sequence[i]->minTau_s;
            g.neutral = true;
            g.sVal = "";
          }
          g.duration_s = sliceDuration_s;
          slice[numSlices][i] = g;
        }

        // Prevent excessive array index
        if (numSlices < MAX_SLICES-1)
        {
          numSlices++;
        }
      }

      pos_s = nearestPos_s;
      gesturePos_s[k] = nearestPos_s;
      index[k]++;
    }
  }




//printf("numSlices=%d\n", numSlices);


  // ****************************************************************
  // Make that neutral consonantal gestures get the value/name of the
  // last previous non-neutral gesture.
  // ****************************************************************

  for (i=1; i <= 3; i++)
  {
    for (k=1; k < numSlices; k++)
    {
      if ((slice[k][i].sVal == "") || (slice[k][i].neutral))
      {
        slice[k][i].sVal = slice[k-1][i].sVal;
        // Make all consonant gestures without a name/value neutral!
        slice[k][i].neutral = true;
      }
    }
  }

  

  // ****************************************************************
  // Calculate the vocal tract parameter targets for each slice.
  // ****************************************************************

  const int NUM_CONSONANT_TYPES = 3;
  double targetValue[VocalTract::NUM_PARAMS];
  double targetTau[VocalTract::NUM_PARAMS];
  
  std::string consonantName[NUM_CONSONANT_TYPES];
  double consonantTau_s[NUM_CONSONANT_TYPES];
  bool isNeutralConsonant[NUM_CONSONANT_TYPES];
  bool consonantExists[NUM_CONSONANT_TYPES];
  double consonantParam[NUM_CONSONANT_TYPES][VocalTract::NUM_PARAMS];
  
  double vowelTau;
  int shapeIndex = 0;
  Target target;

  // Coordinates of the base vowel in the 2D-vowel subspace
  double alphaTongue;
  double betaTongue;
  double alphaLips;
  double betaLips;


  // Clear all targets.
  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    tractParamTargets[i].clear();
  }

  for (i=0; i < numSlices; i++)
  {
    // **************************************************************
    // Set the target values to the vowel gesture shape.
    // **************************************************************

    shapeIndex = vocalTract->getShapeIndex( slice[i][VOWEL_GESTURE].sVal );
    vowelTau = slice[i][VOWEL_GESTURE].tau_s;

    // Try to get the index of Schwa.
    if (shapeIndex == -1)
    {
      shapeIndex = vocalTract->getShapeIndex("@");
    }
    
    if (shapeIndex == -1)
    {
      for (k=0; k < VocalTract::NUM_PARAMS; k++)
      {
        targetValue[k] = vocalTract->param[k].neutral;
        targetTau[k] = vowelTau;
      }
    }
    else
    {
      for (k=0; k < VocalTract::NUM_PARAMS; k++)
      {
        targetValue[k] = vocalTract->shapes[shapeIndex].param[k];
        targetTau[k] = vowelTau;
      }
    }

    // Get the coordinates of the vowel shape in the 2D-vowel subspace.
    mapToVowelSubspace(vocalTract, targetValue, alphaTongue, betaTongue, alphaLips, betaLips);
    limitVowelSubspaceCoord(alphaTongue, betaTongue, alphaLips, betaLips);
    
    // **************************************************************
    // Modify the target values and taus by consonantal gestures.
    // **************************************************************

    consonantName[0] = slice[i][LIP_GESTURE].sVal;
    consonantName[1] = slice[i][TONGUE_TIP_GESTURE].sVal;
    consonantName[2] = slice[i][TONGUE_BODY_GESTURE].sVal;
    consonantTau_s[0] = slice[i][LIP_GESTURE].tau_s;
    consonantTau_s[1] = slice[i][TONGUE_TIP_GESTURE].tau_s;
    consonantTau_s[2] = slice[i][TONGUE_BODY_GESTURE].tau_s;
    isNeutralConsonant[0] = slice[i][LIP_GESTURE].neutral;
    isNeutralConsonant[1] = slice[i][TONGUE_TIP_GESTURE].neutral;
    isNeutralConsonant[2] = slice[i][TONGUE_BODY_GESTURE].neutral;
    consonantExists[0] = false;
    consonantExists[1] = false;
    consonantExists[2] = false;
    
    int numActiveConsonants = 0;
    double consonantTauSum = 0.0;
    double consonantParamSum[VocalTract::NUM_PARAMS] = { 0.0 };

    for (k=0; k < NUM_CONSONANT_TYPES; k++)
    {
      if (isNeutralConsonant[k] == false)
      {
        consonantExists[k] = getContextDependentConsonant(vocalTract, consonantName[k].c_str(), 
          alphaTongue, betaTongue, alphaLips, betaLips, consonantParam[k]);

        // Consonant is not neutral and exists.
        if (consonantExists[k])
        {
          numActiveConsonants++;
          consonantTauSum+= consonantTau_s[k];
          for (m=0; m < VocalTract::NUM_PARAMS; m++)
          {
            consonantParamSum[m]+= consonantParam[k][m];
          }
        }
      }
    }

    // **************************************************************
    // Modify the target values and taus of the vowel when at least 
    // one consonant is active.
    // **************************************************************
    
    if (numActiveConsonants > 0)
    {
      for (m=0; m < VocalTract::NUM_PARAMS; m++)
      {
        targetValue[m] = consonantParamSum[m] / (double)numActiveConsonants;
        targetTau[m] = consonantTauSum / (double)numActiveConsonants;
      }

      // Each consonant has FULL CONTROL of its primary articulator
      
      // Labial consonant:
      if (consonantExists[0])
      {
        targetValue[VocalTract::LP] = consonantParam[0][VocalTract::LP];
        targetValue[VocalTract::LD] = consonantParam[0][VocalTract::LD];
      }

      // Tongue tip consonant:
      if (consonantExists[1])
      {
        targetValue[VocalTract::TTX] = consonantParam[1][VocalTract::TTX];
        targetValue[VocalTract::TTY] = consonantParam[1][VocalTract::TTY];
        targetValue[VocalTract::TS3] = consonantParam[1][VocalTract::TS3];
      }

      // Tongue body consonant:
      if (consonantExists[2])
      {
        targetValue[VocalTract::TCX] = consonantParam[2][VocalTract::TCX];
        targetValue[VocalTract::TCY] = consonantParam[2][VocalTract::TCY];
        targetValue[VocalTract::TS2] = consonantParam[2][VocalTract::TS2];
      }
    }


    // **************************************************************
    // Overwrite the velum targets by the velum gestures.
    // **************************************************************

    if (slice[i][VELIC_GESTURE].neutral == false)
    {
      targetValue[VocalTract::VO] = slice[i][VELIC_GESTURE].dVal;
      targetTau[VocalTract::VO]   = slice[i][VELIC_GESTURE].tau_s;
    }

    // **************************************************************
    // Add the same targets to the different target sequences.
    // **************************************************************

    for (k=0; k < VocalTract::NUM_PARAMS; k++)
    {
      target.duration = slice[i][VOWEL_GESTURE].duration_s;
      target.value = targetValue[k];
      target.tau_s = targetTau[k];
      target.slope = 0.0;     // VT parameter targets have no slopes.
      tractParamTargets[k].push_back(target);
    }
  }


  // ****************************************************************
  // Apply the intrinsic, direction-dependent velocities of the
  // parameters (articulators) to the time constants of the targets.
  // ****************************************************************

  Target* prevTarget;
  Target* thisTarget;
  double velocityFactor = 1.0;
  const double EPSILON = 0.001;

  for (k = 0; k < VocalTract::NUM_PARAMS; k++)
  {
    int numTargets = (int)tractParamTargets[k].size();

    for (i = 1; i < numTargets; i++)
    {
      thisTarget = &tractParamTargets[k][i];
      prevTarget = &tractParamTargets[k][i - 1];

      // Depending on whether the parameter value increases or
      // decreases, a different velocity factor is applied.
      if (thisTarget->value > prevTarget->value)
      {
        velocityFactor = vocalTract->anatomy.positiveVelocityFactor[k];
      }
      else
      {
        velocityFactor = vocalTract->anatomy.negativeVelocityFactor[k];
      }

      if (fabs(velocityFactor) < EPSILON)
      {
        velocityFactor = EPSILON;
      }

      // For an INCREASE of the velocity, REDUCE the time constant.
      thisTarget->tau_s /= velocityFactor;
    }
  }

}


// ****************************************************************************
// ****************************************************************************

void GesturalScore::calcGlottisParamTargets()
{
  int i, k;
  int numGlottisParams = (int)glottis->controlParam.size();
  Gesture *gesture = NULL;
  Glottis::Shape *shape = NULL;
  GestureSequence *sequence = NULL;
  int numGestures = 0;
  Target target;

  // ****************************************************************
  // Clear all targets.
  // ****************************************************************

  for (i=0; i < numGlottisParams; i++)
  {
    glottisParamTargets[i].clear();
  }

  // ****************************************************************
  // Convert the glottal shape gestures into targets.
  // ****************************************************************

  const int FIRST_SHAPE_PARAM = 2;    // "0" is F0 and "1" is lung pressure.

  sequence = &gestures[GLOTTAL_SHAPE_GESTURE];
  numGestures = sequence->numGestures();

  for (i=0; i < numGestures; i++)
  {
    gesture = sequence->getGesture(i);
    shape = glottis->getShape(gesture->sVal);
    if (shape != NULL)
    {
      for (k=FIRST_SHAPE_PARAM; k < numGlottisParams; k++)
      {
        target.duration = gesture->duration_s;
        target.tau_s = gesture->tau_s;
        target.value = shape->controlParam[k];
        target.slope = 0.0;
        
        glottisParamTargets[k].push_back(target);
      }
    }
    else
    {
      for (k=FIRST_SHAPE_PARAM; k < numGlottisParams; k++)
      {
        target.duration = gesture->duration_s;
        target.tau_s = gesture->tau_s;
        target.value = glottis->controlParam[k].neutral;
        target.slope = 0.0;
        
        glottisParamTargets[k].push_back(target);
      }
    }
  }

  // ****************************************************************
  // Convert the F0 gestures into targets.
  // ****************************************************************

  sequence = &gestures[F0_GESTURE];
  numGestures = sequence->numGestures();

  for (i=0; i < numGestures; i++)
  {
    gesture = sequence->getGesture(i);
    
    target.duration = gesture->duration_s;
    target.tau_s = gesture->tau_s;
    target.value = gesture->dVal;
    target.slope = gesture->slope;
    
    glottisParamTargets[Glottis::FREQUENCY].push_back(target);
  }

  // ****************************************************************
  // Convert the lung pressure gestures into targets.
  // ****************************************************************

  sequence = &gestures[PRESSURE_GESTURE];
  numGestures = sequence->numGestures();

  for (i=0; i < numGestures; i++)
  {
    gesture = sequence->getGesture(i);
    
    target.duration = gesture->duration_s;
    target.tau_s = gesture->tau_s;
    target.value = gesture->dVal;
    target.slope = gesture->slope;
    
    glottisParamTargets[Glottis::PRESSURE].push_back(target);
  }
}


// ****************************************************************************
/// Calculate the filter response analytically and then sample this curve to
/// obtain the paramCurve[] samples using a 5th-order system.
// ****************************************************************************

void GesturalScore::calcParamCurve(vector<Target> &paramTargets, vector<double> &paramCurve)
{
  const double EPSILON = 0.000000001;
  int i;

  if (paramTargets.size() < 1)
  {
    return;
  }

  int numTargets = (int)paramTargets.size();
  int targetIndex = 0;
  Target *target = &paramTargets[targetIndex];
  double targetPos_s = 0.0;
  double t_s;
  double t, t2, t3, t4;
  double a, a2, a3, a4;
  double c0, c1, c2, c3, c4;
  double f0, f1, f2, f3, f4;

  // The last sample to consider.
  int finalPointIndex = (int)((getScoreDuration_s() + 0.010) * CURVE_SAMPLING_RATE);
  if (finalPointIndex >= MAX_CURVE_SAMPLES)
  {
    finalPointIndex = MAX_CURVE_SAMPLES - 1;
  }
  
  if (fabs(target->tau_s) < EPSILON)
  {
    target->tau_s = EPSILON;
  }
  a = -1.0 / target->tau_s;
  a2 = a*a;
  a3 = a2*a;
  a4 = a3*a;

  // Coefficients for the initial target are all zero -> we follow the target exactly.
  c0 = 0.0;
  c1 = 0.0;
  c2 = 0.0;
  c3 = 0.0;
  c4 = 0.0;

  // ****************************************************************
  
  for (i = 0; i <= finalPointIndex; i++)
  {
    t_s = (double)i / CURVE_SAMPLING_RATE;

    while ((t_s > targetPos_s + target->duration) && (targetIndex < numTargets - 1))
    {
      // Calculate y(t) and its derivatives at the end of the prev. target.
      
      t = target->duration;
      t2 = t*t;
      t3 = t2*t;
      t4 = t3*t;

      // It's important to consider the slope for f1!
      f0 = exp(a*t)*((c0)+(c1)*t + (c2)*t2 + (c3)*t3 + (c4)*t4) + (target->value + target->duration*target->slope);
      f1 = exp(a*t)*((c0*a + c1) + (c1*a + 2 * c2)*t + (c2*a + 3 * c3)*t2 + (c3*a + 4 * c4)*t3 + (c4*a)*t4) + target->slope;
      f2 = exp(a*t)*((c0*a2 + 2 * c1*a + 2 * c2) + (c1*a2 + 4 * c2*a + 6 * c3)*t + (c2*a2 + 6 * c3*a + 12 * c4)*t2 + (c3*a2 + 8 * c4*a)*t3 + (c4*a2)*t4);
      f3 = exp(a*t)*((c0*a3 + 3 * c1*a2 + 6 * c2*a + 6 * c3) + (c1*a3 + 6 * c2*a2 + 18 * c3*a + 24 * c4)*t + (c2*a3 + 9 * c3*a2 + 36 * c4*a)*t2 + (c3*a3 + 12 * c4*a2)*t3 + (c4*a3)*t4);
      f4 = exp(a*t)*((c0*a4 + 4 * c1*a3 + 12 * c2*a2 + 24 * c3*a + 24 * c4) + (c1*a4 + 8 * c2*a3 + 36 * c3*a2 + 96 * c4*a)*t + (c2*a4 + 12 * c3*a3 + 72 * c4*a2)*t2 + (c3*a4 + 16 * c4*a3)*t3 + (c4*a4)*t4);

      // Go to the next target.
      
      targetPos_s += target->duration;
      targetIndex++;
      target = &paramTargets[targetIndex];

      // Calc. the coefficients for the next target based on the derivatives
      // at the end of the previous target.

      if (fabs(target->tau_s) < EPSILON)
      {
        target->tau_s = EPSILON;
      }
      a = -1.0 / target->tau_s;
      a2 = a*a;
      a3 = a2*a;
      a4 = a3*a;

      c0 = f0 - target->value;
      c1 = (f1 - c0*a - target->slope);   // Slope must be considered here!
      c2 = (f2 - c0*a2 - c1*a * 2) / 2;
      c3 = (f3 - c0*a3 - c1*a2 * 3 - c2*a * 6) / 6;
      c4 = (f4 - c0*a4 - c1*a3 * 4 - c2*a2 * 12 - c3*a * 24) / 24;
    }

    // Calculate the actual curve value here.

    t = t_s - targetPos_s;    // Time relative to the beginning of the target.
    t2 = t*t;
    t3 = t2*t;
    t4 = t3*t;

    paramCurve[i] = exp(a*t)*((c0)+(c1)*t + (c2)*t2 + (c3)*t3 + (c4)*t4) + (target->value + t*target->slope);
  }
}


// ****************************************************************************
