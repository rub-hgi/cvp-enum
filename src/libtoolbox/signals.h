// ----------------------------------------------------------------------------
// Title      : Signals
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : signals.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     extern declaration for signal handling - used for gracefully shutdown
//     during enumeration, if SIGKILL is send to the process.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#ifndef __SIGNALS_H__
#define __SIGNALS_H__

extern bool got_sigterm;

#endif // __SIGNALS_H__
