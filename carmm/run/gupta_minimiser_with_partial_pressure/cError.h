#ifndef __CERROR_H__
#define __CERROR_H__

#include <string>
#include <iostream>
class cError
{
	std::string m_errorText;
public:
	cError(const char* errorText)
	{
		std::cout << "***\n [ERROR] cError Thrown! text: " << errorText << "\n***\n";
		m_errorText = errorText;
	}

	const char* GetText()
	{
		return m_errorText.c_str();
	}
};
#endif
