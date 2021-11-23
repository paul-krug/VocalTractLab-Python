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

#ifndef __SEGMENT_SEQUENCE_H__
#define __SEGMENT_SEQUENCE_H__

#include <string>
#include <vector>

using namespace std;


// ****************************************************************************
/// This class represents a list of key-value pairs (labels) for a phonetic
/// segment.
// ****************************************************************************

class Segment
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  static const int NAME_INDEX = 0;
  static const int DURATION_INDEX = 1;
  static const int MAX_LABELS = 256;

  double duration_s;
  string key[MAX_LABELS];
  string value[MAX_LABELS];
  static const string fixedKey[MAX_LABELS];

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  Segment();
  void reset();
  string getValue(const string &keyName);
  bool setValue(const string &keyName, const string &valueName);
  bool parse(const string &textLine);
  string getTextLine();
  string check();
};


// ****************************************************************************
/// This classes represents a sequence of phonetic segments with several
/// key-value pairs (labels) per segment.
/// This allows the segmentation and labeling of acoustic data.
// ****************************************************************************

class SegmentSequence
{
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  SegmentSequence();
  void clear();
  int numSegments();
  Segment *getSegment(int index);
  bool isValidIndex(int index);

  double getSegmentBegin_s(int index);
  double getSegmentEnd_s(int index);
  int getIndexAt(double pos_s);
  double getDuration_s();
  
  void appendSegment(Segment &s);
  void appendSegment(string name, double duration_s);
  void insertSegment(Segment &s, int index);
  void deleteSegment(int index);
  void setMinSegmentDuration(double duration_s);

  bool writeToFile(const string &fileName);
  bool readFromFile(const string &fileName);
  bool importFromEsps(const string &fileName);

  void resetIteration();
  Segment *getNextPhone(double &startTime_s, double &endTime_s);
  Segment *getNextSyllable(double &startTime_s, double &endTime_s);
  Segment *getNextWord(double &startTime_s, double &endTime_s);
  Segment *getNextPhrase(double &startTime_s, double &endTime_s);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  vector<Segment> segment;
  // Index of the next considered segment for the iteration of units.
  int nextSegment;
  double nextSegmentStartTime_s;
};

#endif