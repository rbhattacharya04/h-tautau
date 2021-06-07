/*! Apply jet uncertainties to the event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <string>
#include <vector>
#include <exception>
#include <TFormula.h>
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "JetCorrectorParameters.h"   // CondFormats/JetMETObjects/interface
#include "JetCorrectionUncertainty.h" // CondFormats/JetMETObjects/interface


enum class Variation {
    NOMINAL = 0,
    DOWN = 1,
    UP = 2
};


template <typename T>
T clip(const T& n, const T& lower, const T& upper) {
    return std::max(lower, std::min(n, upper));
}

namespace JME {

    template <typename T, typename U>
    struct bimap {
        typedef std::unordered_map<T, U> left_type;
        typedef std::unordered_map<U, T> right_type;

        left_type left;
        right_type right;

        bimap(std::initializer_list<typename left_type::value_type> l) {
            for (auto& v: l) {
                left.insert(v);
                right.insert(typename right_type::value_type(v.second, v.first));
            }
        }

        bimap() {
            // Empty
        }

        bimap(bimap&& rhs) {
                    left = std::move(rhs.left);
                    right = std::move(rhs.right);
        }
    };

    enum class Binning {
            JetPt = 0,
            JetEta,
            JetAbsEta,
            JetE,
            JetArea,
            Mu,
            Rho,
            NPV,
    };

};

// Hash function for Binning enum class
namespace std {
    template<>
    struct hash<JME::Binning> {
        typedef JME::Binning argument_type;
        typedef std::size_t result_type;

        hash<uint8_t> int_hash;

        result_type operator()(argument_type const& s) const {
            return int_hash(static_cast<uint8_t>(s));
        }
    };
};

namespace JME {

    class JetParameters {
        public:
            typedef std::unordered_map<Binning, float> value_type;

            JetParameters() = default;
            JetParameters(JetParameters&& rhs);
            JetParameters(std::initializer_list<typename value_type::value_type> init);


            JetParameters& setJetPt(float pt);
            JetParameters& setJetEta(float eta);
            JetParameters& setJetE(float e);
            JetParameters& setJetArea(float area);
            JetParameters& setMu(float mu);
            JetParameters& setRho(float rho);
            JetParameters& setNPV(float npv);
            JetParameters& set(const Binning& bin, float value);
            JetParameters& set(const typename value_type::value_type& value);

            static const bimap<Binning, std::string> binning_to_string;

            std::vector<float> createVector(const std::vector<Binning>& binning) const;

        private:
            value_type m_values;
    };

    class JetResolutionObject {
        public:

            //should I use the std::range or create a range struct just like in the cmssw version

            struct Range {
                float min;
                float max;

                Range() {
                    // Empty
                }

                Range(float min, float max) {
                    this->min = min;
                    this->max = max;
                }

                bool is_inside(float value) const {
                    return (value >= min) && (value <= max);
                }


            };



            class Definition {
                public:
                    Definition() {
                        // Empty
                    }

                    Definition(const std::string& definition);

                    const std::vector<std::string>& getBinsName() const {
                        return m_bins_name;
                    }

                    const std::vector<Binning>& getBins() const {
                        return m_bins;
                    }

                    std::string getBinName(size_t bin) const {
                        return m_bins_name[bin];
                    }

                    size_t nBins() const {
                        return m_bins_name.size();
                    }

                    const std::vector<std::string>& getVariablesName() const {
                        return m_variables_name;
                    }

                    const std::vector<Binning>& getVariables() const {
                        return m_variables;
                    }

                    std::string getVariableName(size_t variable) const {
                        return m_variables_name[variable];
                    }

                    size_t nVariables() const {
                        return m_variables.size();
                    }

                    std::string getFormulaString() const {
                        return m_formula_str;
                    }

                    TFormula* getFormula() const {
                        return m_formula.get();
                    }

                    void init();

                private:

                    std::vector<std::string> m_bins_name;
                    std::vector<std::string> m_variables_name;
                    std::string m_formula_str;

                    std::shared_ptr<TFormula> m_formula;
                    std::vector<Binning> m_bins;
                    std::vector<Binning> m_variables;

            };

            class Record {
                public:

                    Record() {
                        // Empty
                    }

                    Record(const std::string& record, const Definition& def);

                    const std::vector<Range>& getBinsRange() const {
                        return m_bins_range;
                    }

                    const std::vector<Range>& getVariablesRange() const {
                        return m_variables_range;
                    }

                    const std::vector<float>& getParametersValues() const {
                        return m_parameters_values;
                    }

                    size_t nVariables() const {
                        return m_variables_range.size();
                    }

                    size_t nParameters() const {
                        return m_parameters_values.size();
                    }

                private:

                    std::vector<Range> m_bins_range;
                    std::vector<Range> m_variables_range;
                    std::vector<float> m_parameters_values;

            };

            public:

                JetResolutionObject(const std::string& filename);
                //JetResolutionObject(const JetResolutionObject& filename);
                JetResolutionObject();

                void dump() const;
                void saveToFile(const std::string& file) const;

                const Record* getRecord(const JetParameters& bins) const;
                float evaluateFormula(const Record& record, const JetParameters& variables) const;

                const std::vector<Record>& getRecords() const {
                    return m_records;
                }

                const Definition& getDefinition() const {
                    return m_definition;
                }


            private:
                Definition m_definition;
                std::vector<Record> m_records;

                bool m_valid = false;
    };

    class JetResolution {
        public:
            JetResolution(const std::string& filename);

            JetResolution() {
                // Empty
            }

            float getResolution(const JetParameters& parameters) const;

            void dump() const {
                m_object->dump();
            }

            // Advanced usage
            const JetResolutionObject* getResolutionObject() const {
                return m_object.get();
            }

        private:
            std::shared_ptr<JetResolutionObject> m_object;
    };

    class JetResolutionScaleFactor {
        public:
            JetResolutionScaleFactor(const std::string& filename);

            JetResolutionScaleFactor() {
                // Empty
            }

            float getScaleFactor(const JetParameters& parameters, Variation variation = Variation::NOMINAL) const;


            void dump() const {
                m_object->dump();
            }

            // Advanced usage
            const JetResolutionObject* getResolutionObject() const {
                return m_object.get();
            }

        private:
            std::shared_ptr<JetResolutionObject> m_object;

    };

};
