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

#ifndef __SAMPA_H__
#define __SAMPA_H__

#include <string>

using namespace std;

// ****************************************************************************
// ****************************************************************************

class Sampa
{
public:
  static const int NUM_VOWELS = 23;
  static const int NUM_DIPHTHONGS = 24;
  static const int NUM_CONSONANTS = 29;
  static const int NUM_PHONEMES = NUM_VOWELS + NUM_DIPHTHONGS + NUM_CONSONANTS;

  static const string PHONEME[NUM_PHONEMES];

  static int getIndex(const string &symbol);
  static bool isPhoneme(const string &symbol);
  static bool isVowel(const string &symbol);
  static bool isDiphthong(const string& symbol);
  static bool isConsonant(const string &symbol);

  static bool isFricative(const string &symbol);
  static bool isPlosive(const string &symbol);
  static bool isNasal(const string &symbol);
  static bool isLateral(const string &symbol);

  static bool isAlveolar(const string& symbol);
  static bool isVoiced(const string& symbol);
};

#endif
