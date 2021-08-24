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

#include "XmlNode.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>


// ****************************************************************************
// Introduce some static functions.
// ****************************************************************************

static XmlNode *parseElement(const string &input, int &pos, XmlNode *parent, vector<XmlError> *errors);
static void parseStartTag(const string &input, int startPos, int endPos, XmlNode *element, vector<XmlError> *errors);
static string getNextToken(const string &input, int &pos, int endPos, vector<XmlError> *errors);
static bool isValidName(const string &st);
static bool isWhiteSpaceChar(char ch);
static string condenseWhiteSpace(const string &input);
static void decodeErrorPositions(const string &input, vector<XmlError> &errors);


// ****************************************************************************
/// Parse the given XML input text, create an XML-tree, and return the root
/// element of the tree. The first element with the given tag is returned
/// as root node.
/// If an error occured or the input text has no root element with the given
/// tag name, NULL is returned.
/// If a valid XmlNode is returned, it must be deleted by the caller to free
/// the memory.
// ****************************************************************************

XmlNode *xmlParseString(const string &input, const string &tag, vector<XmlError> *errors)
{
  int pos = 0;
  XmlNode *node = NULL;
  bool checkNextElement = true;

  while (checkNextElement)
  {
    node = parseElement(input, pos, NULL, errors);
    if (node == NULL)
    {
      checkNextElement = false;
    }
    else
    {
      if (node->name == tag)
      {
        checkNextElement = false;
      }
      else
      {
        // This is not the element tag that we want!
        delete node;
        node = NULL;
      }
    }
  }

  // ****************************************************************
  // If a root element with the given name was not found -> error.
  // ****************************************************************

  if ((node == NULL) && (errors != NULL))
  {
    XmlError e;
    e.line = 0;
    e.column = 0;
    e.text = "No valid root element '" + tag + "' found!";
    errors->push_back(e);
  }

  // ****************************************************************
  // Calculate the line and column numbers of the error messages
  // ****************************************************************

  if (errors != NULL)
  {
    decodeErrorPositions(input, *errors);
    
    // In the case of any errors, always return NULL.
    if ((errors->size() > 0) && (node != NULL))
    {
      delete node;
      node = NULL;
    }
  }

  return node;
}


// ****************************************************************************
/// Same as xmlParseString(...), except that a file instead of a string is
/// parsed.
/// If a valid XmlNode is returned, it must be deleted by the caller to free
/// the memory.
// ****************************************************************************

XmlNode *xmlParseFile(const string &fileName, const string &tag, vector<XmlError> *errors)
{
  ifstream file(fileName.c_str());
  if (file)
  {
    string text;
    string buffer;

    while (!file.eof())
    {
      getline(file, buffer);
      text+= buffer + "\n";
    }
    file.close();

    return xmlParseString(text, tag, errors);
  }
  else
  {
    printf("Error: File %s could not be opened!\n", fileName.c_str());
    return NULL;
  }
}


// ****************************************************************************
/// Prints the list with errors encountered during parsing an XML document.
// ****************************************************************************

void xmlPrintErrors(vector<XmlError> &errors)
{
  int i;

  if (errors.size() > 0)
  {
    printf("====== XML errors ======\n");

    for (i=0; i < (int)errors.size(); i++)
    {
      printf("ln %d,col %d: %s\n", errors[i].line, errors[i].column, errors[i].text.c_str() );
    }
    printf("\n");
  }
}


// ****************************************************************************
/// Unit test.
// ****************************************************************************

void xmlTest()
{
  vector<XmlError> errors;

  string input =
    "<glottis type=\"Titze\" selected=\"1\">\n"
    "  <static-params>\n"
    "    <param index=\"0\" name=\"Cord rest thickness\" /> \n"
    "    <param index=\"1\" name=\"Cord rest length\" /> \n"
    "  </static-params>\n"
    "  <control-params>\n"
    "    <param index=\"0\" name=\"f0\" /> \n"
    "    <param index=\"1\" name=\"Subglottal pressure\" /> \n"
    "  </control-params>\n"
    "  <shapes>\n"
    "    <shape name=\"default\">\n"
    "      <control-param index=\"0\" value=\"120.000000\" /> \n"
    "      <control-param index=\"1\" value=\"800.000000\" /> \n"
    "    </shape>\n"
    "  </shapes>\n"
    "</glottis>\n";

  printf("Original XML-string\n");
  printf("===================\n\n");

  printf("%s", input.c_str());

  // ****************************************************************
  // Parse and reconstruct the original XML input text.
  // ****************************************************************

  XmlNode *rootNode = xmlParseString(input, "glottis", &errors);

  printf("\nReconstructed XML-string\n");
  printf("========================\n\n");


  if (rootNode != NULL)
  {
    printf("%s", rootNode->toXmlString().c_str());
  }
  else
  {
    printf("No root node returned.\n");
  }

  // ****************************************************************
  // Print the list with errors.
  // ****************************************************************

  if (errors.size() > 0)
  {
    xmlPrintErrors(errors);
  }
}


// ****************************************************************************
/// This function must be called when pos points at a position outside a
/// <...>-tag. When the function is called, pos must point at or before the
/// left angle bracket of the start tag (<tagName ...>).
/// When the function returns, pos points one character after the right angle
/// bracket of the end tag (</tagName>) or the empty tag <tagName .../>.
// ****************************************************************************

static XmlNode *parseElement(const string &input, int &pos, XmlNode *parent, vector<XmlError> *errors)
{
  int i;
  int inputLength = (int)input.length();

  // ****************************************************************
  // Find the opening left angle bracket of the start tag.
  // ****************************************************************

  while ((pos < inputLength) && (input[pos] != '<'))
  {
    pos++;
  }

  if (pos >= inputLength)
  {
    return NULL;
  }

  // Create a new node for the new element
  XmlNode *node = new XmlNode( XmlNode::ELEMENT, parent );

  int startTagLeftPos = pos;
  pos++;              // Move just behind the opening angle bracket

  // Check if this is on of the special elements that don't have an
  // end-tag, like
  // <![CDATA[ ... ]]> or
  // <!DOCTYPE ... > or
  // <?xml ... > or
  // <!-- ... >

  bool hasEndTag = true;
  if ((pos < inputLength) && ((input[pos] == '!') || (input[pos] == '?')))
  {
    hasEndTag = false;
    node->type = XmlNode::OTHER;
  }

  // ****************************************************************
  // Find the closing right angle bracket of the start tag.
  // ****************************************************************

  while ((pos < inputLength) && (input[pos] != '>'))
  {
    // There is a nested node in the start-tag.
    if (input[pos] == '<')
    {
      XmlNode *childNode = parseElement(input, pos, node, errors);
      if (childNode != NULL)
      {
        node->child.push_back(childNode);
      }
    }
    else
    {
      pos++;
    }
  }

  // Error: No right angle bracket closing the start-tag found!
  if (pos >= inputLength)
  {
    if (errors != NULL)
    {
      XmlError e;
      e.line = startTagLeftPos;
      e.column = 0;
      e.text = "Closing bracket '>' missing for start-tag!";
      errors->push_back(e);
    }
    delete node;
    return NULL;
  }

  int startTagRightPos = pos;

  // ****************************************************************
  // Parse the start tag for the name and the attributes of the 
  // element, if it is a regular element.
  // If it is a special element (type=OTHER), put the characters of
  // the element in the text field.
  // ****************************************************************

  if (node->type == XmlNode::ELEMENT)
  {
    parseStartTag(input, startTagLeftPos, startTagRightPos, node, errors);
  }
  else
  if (node->type == XmlNode::OTHER)
  {
    // Names and attributes of the node stay empty.
    node->text = input.substr(startTagLeftPos, startTagRightPos - startTagLeftPos + 1);
  }

  // ****************************************************************
  // Is the tag already closed at the end with "/>" ?
  // If so, we must not search for an end-tag and can return the node.
  // ****************************************************************

  if ((pos > 0) && (input[pos-1] == '/'))
  {
    hasEndTag = false;
  }

  // Go to one position behind the closing bracket of the start-tag
  // before we possibly leave the function.
  pos++;

  // This must be a separate check!
  if (hasEndTag == false)
  {
    return node;
  }

  // ****************************************************************
  // This element should have an end tag. Let's search for it.
  // Thereby, collect the content text(s) and sub-elements of the 
  // element.
  // ****************************************************************

  string text;
  text.reserve(1024);   // Reserve some space to save time.

  while ((pos+1 < inputLength) && 
         ((input[pos] != '<') || (input[pos+1] != '/')))
  {
    // There is a child element beginning in the content part.

    if (input[pos] == '<')
    {
      // Check if we have to write out a text child node.
      text = condenseWhiteSpace(text);
      if (text.empty() == false)
      {
        XmlNode *textNode = new XmlNode( XmlNode::TEXT, node );
        textNode->text = text;
        node->child.push_back(textNode);
        text.clear();
      }

      // Parse the child element node.
      XmlNode *childNode = parseElement(input, pos, node, errors);
      if (childNode != NULL)
      {
        node->child.push_back(childNode);
      }
    }
    else

    // Collect normal text chars for text child nodes.
    {
      text+= input[pos];
      pos++;
    }
  }

  // Check if we have to write out a final text child node.

  text = condenseWhiteSpace(text);
  if (text.empty() == false)
  {
    XmlNode *textNode = new XmlNode( XmlNode::TEXT, node );
    textNode->text = text;
    node->child.push_back(textNode);
  }

  // Make a separate list of childs that contain only sub-elements

  for (i=0; i < (int)node->child.size(); i++)
  {
    if (node->child[i]->type == XmlNode::ELEMENT)
    {
      node->childElement.push_back( node->child[i] );
    }
  }


  // Error: No end-tag of the element found.
  if (pos+1 >= inputLength)
  {
    if (errors != NULL)
    {
      XmlError e;
      e.line = startTagLeftPos;
      e.column = 0;
      e.text = "End-tag '</...>' missing for the element '" + node->name + "'!";
      errors->push_back(e);
    }
    delete node;
    return NULL;
  }

  int endTagLeftPos = pos;
  pos++;

  // ****************************************************************
  // Find the closing right angle bracket of the end tag.
  // ****************************************************************

  while ((pos < inputLength) && (input[pos] != '>'))
  {
    pos++;
  }

  // Error: No right angle bracket closing the end-tag found!
  if (pos >= inputLength)
  {
    if (errors != NULL)
    {
      XmlError e;
      e.line = endTagLeftPos;
      e.column = 0;
      e.text = "Closing bracket '>' missing for end-tag!";
      errors->push_back(e);
    }
    delete node;
    return NULL;
  }

  int endTagRightPos = pos;
  // Go to the position right after the closing bracket
  pos++;


  // ****************************************************************
  // Check if the names in the start tag and the end tag are equal!
  // ****************************************************************

  // Skip the '</' sequence at the beginning of the end tag
  int x = endTagLeftPos + 2;    
  string endTagName = getNextToken(input, x, endTagRightPos - 1, errors);

  if (endTagName != node->name)
  {
    if (errors != NULL)
    {
      XmlError e;
      e.line = endTagLeftPos;
      e.column = 0;
      e.text = "Element names are not equal in the start-tag and the end-tag ('" + 
        node->name + "' vs. '" + endTagName + "') !";
      errors->push_back(e);
    }
    delete node;
    return NULL;
  }

  // Everything ok -> return the created node.
  return node;
}


// ****************************************************************************
/// Parse the given start tag in the input string and set the name and the 
/// attributes of the given element correspondingly.
// ****************************************************************************

static void parseStartTag(const string &input, int startPos, int endPos, XmlNode *element, vector<XmlError> *errors)
{
  // Skip the leading '<' character.
  int pos = startPos + 1;

  // ****************************************************************
  // Parse the name of the element.
  // ****************************************************************

  element->name = getNextToken(input, pos, endPos - 1, errors);
  
  if (isValidName(element->name) == false)
  {
    if (errors != NULL)
    {
      XmlError e;
      e.line = startPos + 1;
      e.column = 0;
      e.text = "'" + element->name + "' is not a valid element name!";
      errors->push_back(e);
    }
  }

  // ****************************************************************
  // Parse the attribute-value pairs.
  // ****************************************************************

  string attName;
  string attValue;
  string equalSign;

  while ((attName = getNextToken(input, pos, endPos - 1, errors)).empty() == false)
  {
    if ((isValidName(attName) == false) && (attName != "/"))
    {
      if (errors != NULL)
      {
        XmlError e;
        e.line = pos;
        e.column = 0;
        e.text = "'" + attName + "' is not a valid attribute name!";
        errors->push_back(e);
      }
      break;
    }

    equalSign = getNextToken(input, pos, endPos - 1, errors);
    attValue  = getNextToken(input, pos, endPos - 1, errors);

    if (equalSign == "=")
    {
      XmlAttribute att;
      att.name = attName;
      att.value = attValue;
      element->attribute.push_back(att);
    }
    else

    // In empty elements, '/' just before the closing '>' is errorly
    // detected as attribute name.
    if (attName != "/")
    {
      if (errors != NULL)
      {
        XmlError e;
        e.line = pos;
        e.column = 0;
        e.text = "Invalid attribute definition (" + attName + " : " + equalSign + " : " + attValue + ")!";
        errors->push_back(e);
      }
      break;
    }
  }
}


// ****************************************************************************
/// This function returns the next token between the given position and
/// the end position in a start-tag or end-tag.
/// Tokens are:
/// o A regular element name.
/// o A string enclosed by single quotes ('...'). 
/// o A string enclosed by double quotes ("..."). 
/// o Special single characters (= or < or > or /).
///
/// White spaces between tokens are ignored (but not white spaces within quotes). 
/// Tokens within quotes in the input string are returned without the quotes.
/// The position pos is moved right behind the returned token.
// ****************************************************************************

static string getNextToken(const string &input, int &pos, int endPos, vector<XmlError> *errors)
{
  // Move to the next non-white-space char.
  while ((pos <= endPos) && (isWhiteSpaceChar(input[pos])))
  {
    pos++;
  }

  if (pos > endPos)
  {
    return string("");
  }

  string token;
  // Reserve some space for the expected characters for performance reasons
  token.reserve(256);   

  // ****************************************************************
  // The token is the equality sign.
  // ****************************************************************

  if ((input[pos] == '=') || (input[pos] == '/') || 
      (input[pos] == '<') || (input[pos] == '>'))
  {
    token = input[pos];
    pos++;
  }
  else

  // ****************************************************************
  // The token is between single quotes.
  // ****************************************************************

  if (input[pos] == '\'')
  {
    // Got to the char after the leading quote.
    pos++;
    // Collect all characters up to the next single quote.
    while ((pos <= endPos) && (input[pos] != '\''))
    {
      token+= input[pos];
      pos++;
    }

    if ((pos > endPos) && (errors != NULL))
    {
      XmlError e;
      e.line = pos;
      e.column = 0;
      e.text = "End quote (') for attribute value missing!";
      errors->push_back(e);
    }

    // Got to the char after the trailing quote.
    pos++;
  }
  else

  // ****************************************************************
  // The token is between double quotes.
  // ****************************************************************

  if (input[pos] == '"')
  {
    // Got to the char after the leading quote.
    pos++;
    // Collect all characters up to the next qouble quote.
    while ((pos <= endPos) && (input[pos] != '"'))
    {
      token+= input[pos];
      pos++;
    }

    if ((pos > endPos) && (errors != NULL))
    {
      XmlError e;
      e.line = pos;
      e.column = 0;
      e.text = "End quote (\") for attribute value missing!";
      errors->push_back(e);
    }

    // Got to the char after the trailing quote.
    pos++;
  }
  else

  // ****************************************************************
  // The token is a regular name.
  // ****************************************************************
  {
    // Collect chars until we find the next white space or reserved char
    while ((pos <= endPos) && 
      (isWhiteSpaceChar(input[pos]) == false) && 
      (input[pos] != '=') &&
      (input[pos] != '<') &&
      (input[pos] != '>') &&
      (input[pos] != '/'))
    {
      token+= input[pos];
      pos++;
    }
  }

  return token;
}


// ****************************************************************************
/// Returns true, if the given string is a valid element name.
// ****************************************************************************

static bool isValidName(const string &st)
{
  int length = (int)st.length();
  if (length < 1)
  {
    return false;
  }

  // ****************************************************************
  // Check the start character.
  // ****************************************************************

  char ch = st[0];
  if ((ch == ':') ||
      (ch == '_') ||
      ((ch >= 'A') && (ch <= 'Z')) ||
      ((ch >= 'a') && (ch <= 'z')))
  {
    // Char is ok.
  }
  else
  {
    return false;
  }

  // ****************************************************************
  // Check the other characters.
  // ****************************************************************

  int i;

  for (i=0; i < length; i++)
  {
    if ((ch == ':') ||
        (ch == '_') ||
        (ch == '-') ||
        (ch == '.') ||
        ((ch >= '0') && (ch <= '9')) ||
        ((ch >= 'A') && (ch <= 'Z')) ||
        ((ch >= 'a') && (ch <= 'z')))
    {
      // Char is ok.
    }
    else
    {
      return false;
    }
  }

  return true;
}

// ****************************************************************************
/// Returns true, if the given character ch is a white space character
/// (either space, carriage return, line feed, tab).
// ****************************************************************************

static bool isWhiteSpaceChar(char ch)
{
  if ((ch == ' ') || (ch == '\x09') || (ch == '\x0d') || (ch == '\x0a'))
  {
    return true;
  }
  else
  {
    return false;
  }
}


// ****************************************************************************
/// This function reduces all sequences of two or more white space characters
/// between normal text characters to exactly one space character.
/// White space at the beginning and the end of the input text is completely
/// removed.
// ****************************************************************************

static string condenseWhiteSpace(const string &input)
{
  int i;
  int inputLength = (int)input.length();
  bool normalCharSeen = false;
  
  string output;
  output.reserve(1024);

  for (i=0; i < inputLength; i++)
  {
    // We are on a non-white-space char.
    if (isWhiteSpaceChar( input[i] ) == false)
    {
      if (i == 0)
      {
        output+= input[i];
      }
      else
      {
        // Add a single space char before the text char.
        if ((isWhiteSpaceChar( input[i-1] )) && (normalCharSeen))
        {
          output+= ' ';
        }
        output+= input[i];
      }

      normalCharSeen = true;
    }
  }

  return output;
}


// ****************************************************************************
/// This function transforms the linear character positions of the errors into
/// line and column numbers for the given input string.
// ****************************************************************************

static void decodeErrorPositions(const string &input, vector<XmlError> &errors)
{
  int i;

  // ****************************************************************
  // Put the linear char positions of new line starts in a list.
  // ****************************************************************

  vector<int> lineStartPos;
  lineStartPos.reserve(1024);   // Maximal expected number of lines
  lineStartPos.push_back(0);    // First line starts at absolute position 0.

  int inputLength = (int)input.length();

  for (i=0; i < inputLength; i++)
  {
    if (input[i] == '\n')
    {
      if ((i < inputLength-1) && (input[i+1] == '\r'))
      {
        lineStartPos.push_back(i+2);
      }
      else
      {
        lineStartPos.push_back(i+1);
      }
    }
  }

  // ****************************************************************
  // Transform the linear absolute positions in line and column 
  // numbers.
  // ****************************************************************

  int numLines = (int)lineStartPos.size();
  int pos, line;

  for (i=0; i < (int)errors.size(); i++)
  {
    pos = errors[i].line;

    // Find the line
    line = 0;
    while ((line < numLines) && (lineStartPos[line] <= pos))
    {
      line++;
    }
    line--;   // Go back one line

    if (line < 0)
    {
      line = 0;
    }
    if (line >= numLines)
    {
      line = numLines - 1;
    }

    // Start counting lines and columns with "1" instead of "0".
    errors[i].line = line + 1;      
    errors[i].column = pos - lineStartPos[line] + 1;
  }
}


// ****************************************************************************
/// Returns the string containing the XML-text corresponding to this element
/// and all of its childs.
// ****************************************************************************

string XmlNode::toXmlString()
{
  ostringstream os;
  toXmlString(os, 0);

  return os.str();
}


// ****************************************************************************
/// Constructor. Constructs an empty element.
// ****************************************************************************

XmlNode::XmlNode(NodeType t, XmlNode *parent)
{
  type = t;
  this->parent = parent;
}

// ****************************************************************************
/// Clear this node and all of its children.
// ****************************************************************************

XmlNode::~XmlNode()
{
  int i;
  for (i=0; i < (int)child.size(); i++)
  {
    delete child[i];
  }
}


// ****************************************************************************
/// Returns the number of child elements with the given name.
// ****************************************************************************

int XmlNode::numChildElements(const string &name)
{
  int count = 0;
  int i;

  for (i=0; i < (int)childElement.size(); i++)
  {
    if (childElement[i]->name == name)
    {
      count++;
    }
  }

  return count;
}


// ****************************************************************************
/// Returns the child element with the given name and index, or NULL, if it 
/// doesn't exist.
// ****************************************************************************

XmlNode *XmlNode::getChildElement(const string &name, int index)
{
  int counter = 0;
  int i;
  int k = -1;
  
  for (i=0; i < (int)childElement.size(); i++)
  {
    if (childElement[i]->name == name)
    {
      if (index == counter)
      {
        k = i;
        break;
      }
      counter++;
    }
  }

  if (k != -1)
  {
    return childElement[k];
  }
  else
  {
    return NULL;
  }
}


// ****************************************************************************
/// Returns true, if this node has an attribute with the given name, and 
/// otherwise false.
// ****************************************************************************

bool XmlNode::hasAttribute(const string &name)
{
  int i;

  for (i=0; i < (int)attribute.size(); i++)
  {
    if (attribute[i].name == name)
    {
      return true;
    }
  }

  return false;
}


// ****************************************************************************
/// Returns the integer value of the given attribute, or 0, if the attribute
/// does not exist or is not convertable to integer.
// ****************************************************************************

int XmlNode::getAttributeInt(const string &name)
{
  int i;
  int k = -1;

  for (i=0; i < (int)attribute.size(); i++)
  {
    if (attribute[i].name == name)
    {
      k = i;
      break;
    }
  }

  if (k != -1)
  {
    return atoi( attribute[k].value.c_str() );
  }
  else
  {
    return 0;
  }
}


// ****************************************************************************
/// Returns the double value of the given attribute, or 0.0, if the attribute
/// does not exist or is not convertable to double.
// ****************************************************************************

double XmlNode::getAttributeDouble(const string &name)
{
  int i;
  int k = -1;

  for (i=0; i < (int)attribute.size(); i++)
  {
    if (attribute[i].name == name)
    {
      k = i;
      break;
    }
  }

  if (k != -1)
  {
    return atof( attribute[k].value.c_str() );
  }
  else
  {
    return 0.0;
  }
}


// ****************************************************************************
/// Returns the string value of the given attribute, or an empty string, if the
/// attribute does not exist.
// ****************************************************************************

string XmlNode::getAttributeString(const string &name)
{
  int i;
  int k = -1;

  for (i=0; i < (int)attribute.size(); i++)
  {
    if (attribute[i].name == name)
    {
      k = i;
      break;
    }
  }

  if (k != -1)
  {
    return attribute[k].value;
  }
  else
  {
    return string("");
  }
}


// ****************************************************************************
/// This function adds the XML-text corresponding to this element and all its
/// child elements to the output stream.
// ****************************************************************************

void XmlNode::toXmlString(ostream &os, int indent)
{
  char st[1024];
  int i;

  // ****************************************************************
  // This is a regular element.
  // ****************************************************************

  if (type == ELEMENT)
  {
    // Output the start tag of the element.
    os << string(indent, ' ') << "<" << name;

    for (i=0; i < (int)attribute.size(); i++)
    {
      sprintf(st, " %s=\"%s\"", attribute[i].name.c_str(), attribute[i].value.c_str());
      os << st;
    }

    // If the element has no childs and no content text, 
    // close the element and return.
    if ((child.size() == 0) && (text.size() == 0))
    {
      os << "/>" << endl;
      return;
    }
    else
    {
      os << ">" << endl;
    }

    // Output the child elements and the text content.
    for (i=0; i < (int)child.size(); i++)
    {
      child[i]->toXmlString(os, indent + 2);
      
      // Check if parent attributes are set correctly...
      if (child[i]->parent != this)
      {
        os << "ERROR: PARENT MEMBER VARIABLE INCORRECT!" << endl;
      }
    }

    // Output the end tag of the element.
    os << string(indent, ' ') << "</" << name << ">" << endl;
  }
  else

  // ****************************************************************
  // This is a text section (a leaf node).
  // ****************************************************************

  if (type == TEXT)
  {
    os << string(indent, ' ') << text << endl;
  }
  else

  // ****************************************************************
  // This is an other leave node (comment, CDATA, etc.)
  // ****************************************************************

  if (type == OTHER)
  {
    os << string(indent, ' ') << text << endl;
  }
}

// ****************************************************************************

