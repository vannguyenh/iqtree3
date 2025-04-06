/*
 * modelmultistate.h
 *
 *  Created on: Apr 6, 2025
 *      Author: Hiroaki Sato
 *      Original author: minh
 */

#ifndef MODELMULTISTATE_H_
#define MODELMULTISTATE_H_

#include "modelmarkov.h"

class ModelMultistate: public ModelMarkov {
public:
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelMultistate(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree);


	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

    /**
     return the number of dimensions
     */
    virtual int getNDim();

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();
    
    /**
     save object into the checkpoint
     */
    virtual void saveCheckpoint();
    
    /**
     restore object from the checkpoint
     */
    virtual void restoreCheckpoint();

    /**
     * @return model name
     */
    virtual string getName();

    /**
     * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
     */
    virtual string getNameParams(bool show_fixed_params = false);

    /**
     write information
     @param out output stream
     */
    virtual void writeInfo(ostream &out);

    /**
     write parameters, used with modeltest
     @param out output stream
     */
    virtual void writeParameters(ostream &out);

	/**
		read the rates from an input stream. it will throw error messages if failed
		@param in input stream
	*/
    virtual void readRates(istream &in) noexcept(false);

		virtual ~ModelMultistate();
};

#endif /* MODELMUTISTATE_H_ */
