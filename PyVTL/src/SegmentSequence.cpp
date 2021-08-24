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

#include "SegmentSequence.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

// ****************************************************************************
// ****************************************************************************

const string Segment::fixedKey[Segment::MAX_LABELS] =
{ 
  // Phone-level parameters
  "name", 
  "duration_s",
  // Syllable-level parameters
  "start_of_syllable", 
  "word_accent", 
  "phrase_accent",
  "pitch_target_offset_st",
  "pitch_target_slope_st_s",
  // Word-level parameters
  "start_of_word",
  "word_orthographic",
  "word_canonic",
  "part_of_speech",
  // Intonation-phrase-level parameters
  "start_of_phrase",
  // Sentence-level parameters
  "start_of_sentence",
  "sentence_type",
  ""
};


// ****************************************************************************
/// Constructor.
// ****************************************************************************

Segment::Segment()
{
  reset();
}


// ****************************************************************************
/// Resets the segment labels.
// ****************************************************************************

void Segment::reset()
{
  int i;
  for (i=0; i < MAX_LABELS; i++)
  {
    key[i] = fixedKey[i];
    value[i] = "";
  }

  duration_s = 0.0;
}


// ****************************************************************************
/// Returns the value for the given key, or an empty string, if the key does
/// not exist.
// ****************************************************************************

string Segment::getValue(const string &keyName)
{
  int i, k = -1;
  for (i=0; i < MAX_LABELS; i++)
  {
    if (key[i] == keyName)
    {
      k = i;
      break;
    }
  }

  if (k != -1)
  {
    return value[k];
  }
  else
  {
    return string("");
  }
}


// ****************************************************************************
/// Sets the given value for the given key. If the key does not exist, false
/// is returned, and otherwise true.
// ****************************************************************************

bool Segment::setValue(const string &keyName, const string &valueName)
{
  int i, k = -1;
  for (i=0; i < MAX_LABELS; i++)
  {
    if (key[i] == keyName)
    {
      k = i;
      break;
    }
  }

  if (k != -1)
  {
    value[k] = valueName;
    return true;
  }
  else
  {
    printf("Segment::setValue(): The key %s does not exist!\n", keyName.c_str());
    return false;
  }
}


// ****************************************************************************
/// Parse the given text line with respect to the key-value pairs.
/// The syntax of the line should be:
/// key1 = value1; key2 = value2; ...
/// If no (key, value) pairs were read, false is returned.
// ****************************************************************************

bool Segment::parse(const string &line)
{
  string readKey, readValue;
  int pos = 0;
  int length = (int)line.length();
  int i, k;
  int freeLabelIndex;
  int numPairsRead = 0;

  // Clear the current segment data.

  reset();

  // ****************************************************************
  // Retrieve one key-value pair in each loop.
  // ****************************************************************

  while (pos < length)
  {
    // Skip all white spaces.
    while ((pos < length) && ((line[pos] == ' ') || (line[pos] == '\t'))) { pos++; }

    // Read the key (up to the next "=" sign).
    readKey = "";
    while ((pos < length) && (line[pos] != '=')) 
    { 
      readKey+= line[pos];
      pos++; 
    }

    // Got to the next char after the "=".
    pos++;

    // Skip all white spaces.
    while ((pos < length) && ((line[pos] == ' ') || (line[pos] == '\t'))) { pos++; }

    // Read the value (up to the next semicolon or the end of the line).
    readValue = "";
    while ((pos < length) && (line[pos] != ';')) 
    { 
      readValue+= line[pos];
      pos++; 
    }

    // Got to the next char after the ";".
    pos++;

    // **************************************************************
    // Remove all white spaces at the end of the key and the value 
    // string.
    // **************************************************************
    
    while (((int)readKey.length() > 0) && (readKey[readKey.length()-1] == ' '))
    {
      readKey = readKey.substr(0, readKey.length()-1);
    }
    while (((int)readValue.length() > 0) && (readValue[readValue.length()-1] == ' '))
    {
      readValue = readValue.substr(0, readValue.length()-1);
    }

    // **************************************************************
    // Put the key-value pair in the list.
    // **************************************************************

    if ((readKey.empty() == false) && (readValue.empty() == false))
    {
      // Is the key one of the fixed keys?
      k = -1;
      freeLabelIndex = -1;
      for (i=0; (i < MAX_LABELS) && (k == -1); i++)
      {
        if ((fixedKey[i].empty() == false) && (fixedKey[i] == readKey))
        {
          k = i;
        }
        if ((key[i].empty()) && (freeLabelIndex == -1))
        {
          freeLabelIndex = i;
        }
      }

      if (k != -1)
      {
        key[k] = fixedKey[k];
        value[k] = readValue;
      }
      else
      if (freeLabelIndex != -1)
      {
        key[freeLabelIndex] = readKey;
        value[freeLabelIndex] = readValue;
      }

      numPairsRead++;
    }

  }

  // Transform the duration string to double
  duration_s = atof(value[DURATION_INDEX].c_str());

  if (numPairsRead > 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}


// ****************************************************************************
/// Puts all key-value pairs into one line of text.
// ****************************************************************************

string Segment::getTextLine()
{
  int i;
  string line;
  char st[1024];

  sprintf(st, "%f", duration_s);
  value[DURATION_INDEX] = string(st);

  for (i=0; i < MAX_LABELS; i++)
  {
    if (((key[i].empty() == false) && (value[i].empty() == false)) || (key[i] == "name"))
    {
      line = line + key[i] + " = " + value[i] + "; ";
    }
  }
  return line;
}


// ****************************************************************************
/// Performs some checks of plausibility of the key and value strings and 
/// returns an according error string.
/// An empty string is returned when no error was encountered.
// ****************************************************************************

string Segment::check()
{
  int i;
  string st = "";    // Default: no error

  // None of the strings may contain "=" or ";"
  for (i=0; i < MAX_LABELS; i++)
  {
    if (key[i].find("=") != string::npos)
    {
      st = "Error: The key " + key[i] + " contains the character '=' !";
      return st;
    }
    if (key[i].find(";") != string::npos)
    {
      st = "Error: The key " + key[i] + " contains the character ';' !";
      return st;
    }
    
    if (value[i].find("=") != string::npos)
    {
      st = "Error: The value " + value[i] + " contains the character '=' !";
      return st;
    }
    if (value[i].find(";") != string::npos)
    {
      st = "Error: The value " + value[i] + " contains the character ';' !";
      return st;
    }
  }

  // The duration should not be smaller than 1 ms.
  if (duration_s < 0.001)
  {
    st = "Warning: The duration is shorter than 1 ms!";
    return st;
  }

  return string("");
}


// ****************************************************************************
/// Constructor.
// ****************************************************************************

SegmentSequence::SegmentSequence()
{
  clear();
}


// ****************************************************************************
// ****************************************************************************

void SegmentSequence::clear()
{
  segment.clear();
  nextSegment = 0;
}


// ****************************************************************************
/// Returns the number of phones.
// ****************************************************************************

int SegmentSequence::numSegments()
{
  return (int)segment.size();
}


// ****************************************************************************
/// Returns the phone labels for the given phone index or NULL, if the index
/// is out of range.
// ****************************************************************************

Segment *SegmentSequence::getSegment(int index)
{
  if ((index >= 0) && (index < (int)segment.size()))
  {
    return &segment[index];
  }
  else
  {
    return NULL;
  }
}


// ****************************************************************************
/// Returns true, if the given index specifies an existing phone.
// ****************************************************************************

bool SegmentSequence::isValidIndex(int index)
{
  if ((index >= 0) && (index < (int)segment.size()))
  {
    return true;
  }
  else
  {
    return false;
  }
}


// ****************************************************************************
/// Returns the start time of the given segment.
// ****************************************************************************

double SegmentSequence::getSegmentBegin_s(int index)
{
  double t_s = 0.0;
  
  if (isValidIndex(index))
  {
    int i;
    for (i=0; i < index; i++)
    {
      t_s+= segment[i].duration_s;
    }
  }

  return t_s;

}


// ****************************************************************************
/// Returns the end time of the given segment.
// ****************************************************************************

double SegmentSequence::getSegmentEnd_s(int index)
{
  double t_s = 0.0;
  
  if (isValidIndex(index))
  {
    int i;
    for (i=0; i <= index; i++)
    {
      t_s+= segment[i].duration_s;
    }
  }

  return t_s;
}


// ****************************************************************************
/// Returns the index of the segment at the position pos_s or -1, if there is
/// no phone.
// ****************************************************************************

int SegmentSequence::getIndexAt(double pos_s)
{
  int i, k = -1;
  double start_s = 0.0;

  for (i=0; (i < (int)segment.size()) && (k == -1); i++)
  {
    if ((pos_s >= start_s) && (pos_s < start_s + segment[i].duration_s))
    { 
      k = i; 
    }
    start_s+= segment[i].duration_s;
  }

  return k;

}


// ****************************************************************************
/// Returns the duration of this sequence.
// ****************************************************************************

double SegmentSequence::getDuration_s()
{
  int i;
  double duration_s = 0.0;

  for (i=0; i < (int)segment.size(); i++)
  {
    duration_s+= segment[i].duration_s;
  }

  return duration_s;
}


// ****************************************************************************
/// Adds the given segment to the end of the sequence.
// ****************************************************************************

void SegmentSequence::appendSegment(Segment &s)
{
  segment.push_back(s);
}


// ****************************************************************************
/// Adds a segment with the given name and duration to the end of the sequence.
// ****************************************************************************

void SegmentSequence::appendSegment(string name, double duration_s)
{
  Segment s;
  s.value[Segment::NAME_INDEX] = name;
  s.duration_s = duration_s;
  s.value[Segment::DURATION_INDEX] = std::to_string(duration_s);
  segment.push_back(s);
}


// ****************************************************************************
/// Inserts a new segment at the position index into the sequence.
// ****************************************************************************

void SegmentSequence::insertSegment(Segment &s, int index)
{
  if (isValidIndex(index) == false) 
  { 
    return; 
  }
  vector<Segment>::iterator iter = segment.begin();
  advance(iter, index);
  segment.insert(iter, s);
}


// ****************************************************************************
/// Deletes the segment at the position index.
// ****************************************************************************

void SegmentSequence::deleteSegment(int index)
{
  if (isValidIndex(index) == false) 
  { 
    return; 
  }
  vector<Segment>::iterator iter = segment.begin();
  advance(iter, index);
  segment.erase(iter);
}


// ****************************************************************************
/// Make sure that every segment has the provided minDuration_s. The total
/// duration increments that are used to give short segments the minimum 
/// duration are subtracted from longer segment, proportional to their duration.
/// So the total duration of the utterance should not change.
// ****************************************************************************

void SegmentSequence::setMinSegmentDuration(double minDuration_s)
{
  const double DELTA_DURATION_S = 0.005;    // 5 ms
  int i;
  Segment* s = NULL;
  string name;
  int N = (int)segment.size();
  double increment_s = 0.0;
  double decrement_s = 0.0;
  double totalIncrement_s = 0.0;
  double totalDurationOfLongSegments = 0.0;

  // ****************************************************************
  // Ensure the minimum duration for all segments and sum up the
  // durations of the *long* segments.
  // ****************************************************************

  for (i = 0; i < N; i++)
  {
    s = &segment[i];
    name = s->value[Segment::NAME_INDEX];

    // Pause segments and glottal stops are not processed.
    if ((name == "") || (name == "?"))
    {
      continue;   // Continue with next iteration.
    }

    if (s->duration_s < minDuration_s)
    {
      increment_s = minDuration_s - s->duration_s;
      totalIncrement_s += increment_s;
      s->duration_s = minDuration_s;
    }
    else
    if (s->duration_s >= minDuration_s + DELTA_DURATION_S)
    {
      totalDurationOfLongSegments += s->duration_s;
    }
  }

  // ****************************************************************
  // Reduce the length of the *long* segments so that the total 
  // duration of the utterance is the same as before.
  // ****************************************************************

  // Safety check.
  if (totalDurationOfLongSegments < 0.001)
  {
    totalDurationOfLongSegments = 0.001;    // 1 ms
  }

  for (i = 0; i < N; i++)
  {
    s = &segment[i];
    name = s->value[Segment::NAME_INDEX];

    // Pause segments and glottal stops are not processed.
    if ((name == "") || (name == "?"))
    {
      continue;   // Continue with next iteration.
    }

    if (s->duration_s >= minDuration_s + DELTA_DURATION_S)
    {
      decrement_s = totalIncrement_s * (s->duration_s / totalDurationOfLongSegments);
      s->duration_s -= decrement_s;

      // Also the shortened *long* segments may not become shorter
      // than the minimum length.
      if (s->duration_s < minDuration_s)
      {
        s->duration_s = minDuration_s;
      }
    }
  }
}


// ****************************************************************************
/// Write a text file with the segment sequence.
// ****************************************************************************

bool SegmentSequence::writeToFile(const string &fileName)
{
  ofstream os(fileName);
  int i;
  string st;

  if (!os)
  {
    return false;
  }

  for (i=0; i < (int)segment.size(); i++)
  {
    st = segment[i].getTextLine();
    os << st << endl;
  }

  return true;
}


// ****************************************************************************
/// Read a text file with the segment sequence.
// ****************************************************************************

bool SegmentSequence::readFromFile(const string &fileName)
{
  ifstream file(fileName);
  
  if (!file)
  {
    return false;
  }

  char st[2048];
  Segment s;

  segment.clear();
  while (file.eof() == false)
  {
    file.getline(st, 2047);
    if (s.parse(string(st)))
    {
      segment.push_back(s);
    }
  }

  return true;
}


// ****************************************************************************
/// Imports a phoneme seuqence from the ESPS/waves+ file format, that was used
/// in the Berlin Emotion database by Sendlmeier.
// ****************************************************************************

bool SegmentSequence::importFromEsps(const string &fileName)
{
  ifstream file(fileName);
  if (!file)
  {
    printf("Cannot open file %s.\n", fileName.c_str());
    return false;
  }

  // Clear the current segments.
  segment.clear();
 
  enum State { WAIT_FOR_BEGINNING, READ_PHONES };
  const int MAX_TOKENS = 5;
  char line[1024];
  istringstream is;
  string token[MAX_TOKENS];
  Segment s;
  string lastName = "-";
  double pos_s;

  double lastPos_s = 0.0;  
  State state = WAIT_FOR_BEGINNING;

  while (file)
  {
    file.getline(line, 1023);

    is.clear();     // clear eof-flag
    is.str(line);

    is >> token[0];

    if (state == WAIT_FOR_BEGINNING)
    {
      if (token[0] == "#")
      {
        state = READ_PHONES;
      }
    }
    else
    if (state == READ_PHONES)
    {
      is >> token[1];
      is >> token[2];   // Phone name
      is >> token[3];

      if ((token[2].empty() == false) && 
          (token[2][0] != '+') && (token[2][0] != '-'))
      {
        pos_s = atof(token[0].c_str());

        s.reset();
        s.value[Segment::NAME_INDEX] = lastName;
        s.duration_s = pos_s - lastPos_s;
        segment.push_back(s);

        lastPos_s = pos_s;
        lastName = token[2];
      }
    }
  }
  
  file.close();
  return true;
}



// ****************************************************************************
/// Resets the iteration of phones, syllables, words, or phrases.
// ****************************************************************************

void SegmentSequence::resetIteration()
{
  nextSegment = 0;
  nextSegmentStartTime_s = 0.0;
}


// ****************************************************************************
/// Returns the next phone in the sequence, or NULL, if there is no next phone. 
/// To initiate the iteration, call resetIteration().
// ****************************************************************************

Segment *SegmentSequence::getNextPhone(double &startTime_s, double &endTime_s)
{
  startTime_s = 0.0;
  endTime_s = 0.0;
  Segment *s = NULL;

  if ((nextSegment >= 0) && (nextSegment < (int)segment.size()))
  {
    s = &segment[nextSegment];
    startTime_s = nextSegmentStartTime_s;
    endTime_s = startTime_s + s->duration_s;

    nextSegmentStartTime_s+= s->duration_s;
    nextSegment++;
  }

  return s;
}


// ****************************************************************************
/// Returns the next syllable in the sequence in terms of its first segment,
/// start time and duration, or NULL, if there is no next syllable.
/// To initiate the iteration, call resetIteration().
// ****************************************************************************

Segment *SegmentSequence::getNextSyllable(double &startTime_s, double &endTime_s)
{
  const string KEY = "start_of_syllable";
  int numSegments = (int)segment.size();
  startTime_s = 0.0;
  endTime_s = 0.0;
  Segment *s = NULL;

  // Iterate to the next segment where start_of_syllable=1.
  while ((nextSegment < numSegments) && 
         (segment[nextSegment].getValue(KEY) != "1"))
  {
    nextSegmentStartTime_s+= segment[nextSegment].duration_s;
    nextSegment++;
  }

  // Found the beginning of a next syllable?
  if ((nextSegment < numSegments) && (segment[nextSegment].getValue(KEY) == "1"))
  {
    s = &segment[nextSegment];
    startTime_s = nextSegmentStartTime_s;

    // Go one segment further.
    nextSegmentStartTime_s+= segment[nextSegment].duration_s;
    nextSegment++;

    // Iterate to the next segment where start_of_syllable=1 or to 
    // the end of the sequence.
    while ((nextSegment < numSegments) && 
           (segment[nextSegment].getValue(KEY) != "1"))
    {
      nextSegmentStartTime_s+= segment[nextSegment].duration_s;
      nextSegment++;
    }

    endTime_s = nextSegmentStartTime_s;
  }

  return s;
}


// ****************************************************************************
/// Returns the next word in the sequence in terms of its first segment,
/// start time and duration, or NULL, if there is no next word.
/// To initiate the iteration, call resetIteration().
// ****************************************************************************

Segment *SegmentSequence::getNextWord(double &startTime_s, double &endTime_s)
{
  const string KEY = "start_of_word";
  int numSegments = (int)segment.size();
  startTime_s = 0.0;
  endTime_s = 0.0;
  Segment *s = NULL;

  // Iterate to the next segment where start_of_syllable=1.
  while ((nextSegment < numSegments) && 
         (segment[nextSegment].getValue(KEY) != "1"))
  {
    nextSegmentStartTime_s+= segment[nextSegment].duration_s;
    nextSegment++;
  }

  // Found the beginning of a next word?
  if ((nextSegment < numSegments) && (segment[nextSegment].getValue(KEY) == "1"))
  {
    s = &segment[nextSegment];
    startTime_s = nextSegmentStartTime_s;

    // Go one segment further.
    nextSegmentStartTime_s+= segment[nextSegment].duration_s;
    nextSegment++;

    // Iterate to the next segment where start_of_word=1 or to 
    // the end of the sequence.
    while ((nextSegment < numSegments) && 
           (segment[nextSegment].getValue(KEY) != "1"))
    {
      nextSegmentStartTime_s+= segment[nextSegment].duration_s;
      nextSegment++;
    }

    endTime_s = nextSegmentStartTime_s;
  }

  return s;
}


// ****************************************************************************
/// Returns the next phrase in the sequence in terms of its first segment,
/// start time and duration, or NULL, if there is no next phrase.
/// To initiate the iteration, call resetIteration().
// ****************************************************************************

Segment *SegmentSequence::getNextPhrase(double &startTime_s, double &endTime_s)
{
  const string KEY = "start_of_phrase";
  int numSegments = (int)segment.size();
  startTime_s = 0.0;
  endTime_s = 0.0;
  Segment *s = NULL;

  // Iterate to the next segment where start_of_phrase=1.
  while ((nextSegment < numSegments) && 
         (segment[nextSegment].getValue(KEY) != "1"))
  {
    nextSegmentStartTime_s+= segment[nextSegment].duration_s;
    nextSegment++;
  }

  // Found the beginning of a next phrase?
  if ((nextSegment < numSegments) && (segment[nextSegment].getValue(KEY) == "1"))
  {
    s = &segment[nextSegment];
    startTime_s = nextSegmentStartTime_s;

    // Go one segment further.
    nextSegmentStartTime_s+= segment[nextSegment].duration_s;
    nextSegment++;

    // Iterate to the next segment where start_of_phrase=1 or to 
    // the end of the sequence.
    while ((nextSegment < numSegments) && 
           (segment[nextSegment].getValue(KEY) != "1"))
    {
      nextSegmentStartTime_s+= segment[nextSegment].duration_s;
      nextSegment++;
    }

    endTime_s = nextSegmentStartTime_s;
  }

  return s;
}


// ****************************************************************************
