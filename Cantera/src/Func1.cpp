#include "Func1.h"
#include "stringUtils.h"
#include "global.h"

using namespace std;

namespace Cantera {

    static Func1* checkDupl(Func1& f) {
        if (f.parent() != 0) 
            return &f.duplicate();
        else
            return &f;
    }

    Func1& Sin1::derivative() const {
        Func1* c = new Cos1(m_c);
        Func1* r = &newTimesConstFunction(*c, m_c);
#ifdef DEBUG_FUNC
        cout << "Sin1::derivative: \n";
        cout << "function = \'" + write("x") + "\'\n";
        cout << "derivative = \'" + r->write("x") + "\'\n";
#endif 
        return *r;
    } 

    Func1& Cos1::derivative() const {
        Func1* s = new Sin1(m_c);
        Func1* r = &newTimesConstFunction(*s, -m_c);
#ifdef DEBUG_FUNC
        cout << "Cos1::derivative: \n";
        cout << "function = \'" + write("x") + "\'\n";
        cout << "derivative = \'" + r->write("x") + "\'\n";
#endif 
        return *r;
    }

    Func1& Exp1::derivative() const {
        Func1* f = new Exp1(m_c);
        if (m_c != 1.0) 
            return newTimesConstFunction(*f, m_c);
        else
            return *f;
    }

    Func1& Pow1::derivative() const {
        Func1* r;
        if (m_c == 0.0) {
            r = new Const1(0.0);
        }
        else if (m_c == 1.0) {
            r = new Const1(1.0);
        }
        else {
            Func1* f = new Pow1(m_c - 1.0);
            r = &newTimesConstFunction(*f, m_c);
        }
#ifdef DEBUG_FUNC
        cout << "Pow1::derivative: \n";
        cout << "function = \'" + write("x") + "\'\n";
        cout << "derivative = \'" + r->write("x") + "\'\n";
#endif 
        return *r;
    }

    string Func1::write(std::string arg) const {
        return "<unknown " + int2str(ID()) + ">("+arg+")";
    }

    string Sin1::write(string arg) const {
        string c  = "";
        if (m_c != 1.0) c = fp2str(m_c);
        return "\\sin("+c+arg+")";
    }

    string Cos1::write(string arg) const {
        string c  = "";
        if (m_c != 1.0) c = fp2str(m_c);
        return "\\cos("+c+arg+")";
    }

    string Pow1::write(string arg) const {
        //cout << "Pow1" << endl;
        string c  = "";
        if (m_c == 0.5) {
            return "\\sqrt{" + arg + "}";
        } 
        if (m_c == -0.5) {
            return "\\frac{1}{\\sqrt{" + arg + "}}";
        } 
        if (m_c != 1.0) {
            c = fp2str(m_c);
            return "\\left("+arg+"\\right)^{"+c+"}";
        }
        else {
            return arg; 
        }
    }

    string Exp1::write(string arg) const {
        string c  = "";
        if (m_c != 1.0) c = fp2str(m_c);
        return "\\exp("+c+arg+")";
    }

    string Const1::write(string arg) const {
        //cout << "Const1" << endl;
        string c  = "";
        c = fp2str(m_c);
        return c;
    }

    string Ratio1::write(string arg) const {
        //cout << "Ratio1" << endl;
        return "\\frac{" + m_f1->write(arg) + "}{" 
            + m_f2->write(arg) + "}";
    }

    string Product1::write(string arg) const {
        //cout << "Product1" << endl;
        string s = m_f1->write(arg); 
        if (m_f1->order() < order()) s = "\\left(" + s + "\\right)";
        string s2 = m_f2->write(arg); 
        if (m_f2->order() < order()) s2 = "\\left(" + s2 + "\\right)";
        return s + " " + s2;
    }

    string Sum1::write(string arg) const {
        //cout << "Sum1" << endl;
        string s1 = m_f1->write(arg);
        string s2 = m_f2->write(arg);
        if (s2[0] == '-') return s1 + " - " + s2.substr(1,s2.size());
        else return s1 + " + " + s2;
    }

    string Diff1::write(string arg) const {
        //cout << "Diff1" << endl;
        string s1 = m_f1->write(arg);
        string s2 = m_f2->write(arg);
        if (s2[0] == '-') return s1 + " + " + s2.substr(1,s2.size());
        else return s1 + " - " + s2;
    }

    string Composite1::write(string arg) const {
        //cout << "Composite1" << endl;        
        string g = m_f2->write(arg);
        return m_f1->write(g);
    }

    string TimesConstant1::write(string arg) const {
        //cout << "TimesConstant1" << endl;
        string s = m_f1->write(arg);
        if (m_f1->order() < order()) s = "\\left(" + s + "\\right)";
        if (m_c == 1.0) return s;
        if (m_c == -1.0) return "-"+s;
        char n = s[0];
        if (n >= '0' && n <= '9')
            s = "\\left(" + s + "\\right)";
        return fp2str(m_c) + s;
    }

    string PlusConstant1::write(string arg) const {
        //cout << "PlusConstant1" << endl;
        if (m_c == 0.0) return m_f1->write(arg);
        return m_f1->write(arg) + " + " + fp2str(m_c); 
    }

    doublereal Func1::isProportional(TimesConstant1& other) {
        if (isIdentical(other.func1())) return other.c();
        return 0.0;
    }
    doublereal Func1::isProportional(Func1& other) {
        if (isIdentical(other)) return 1.0;
        else return 0.0;
    }

    static bool isConstant(Func1& f) { 
        if (f.ID() == ConstFuncType) 
            return true;
        else
            return false;
    }

    static bool isZero(Func1& f) { 
        if (f.ID() == ConstFuncType && f.c() == 0.0) 
            return true;
        else
            return false;
    }

    static bool isOne(Func1& f) { 
        if (f.ID() == ConstFuncType && f.c() == 1.0) 
            return true;
        else
            return false;
    }

    static bool isTimesConst(Func1& f) { 
        if (f.ID() == TimesConstantFuncType) 
            return true;
        else
            return false;
    }

    static bool isExp(Func1& f) { 
        if (f.ID() == ExpFuncType) 
            return true;
        else
            return false;
    }

    static bool isPow(Func1& f) { 
        if (f.ID() == PowFuncType) 
            return true;
        else
            return false;
    }

    Func1& newSumFunction(Func1& f1, Func1& f2) {
        if (f1.isIdentical(f2)) 
            return newTimesConstFunction(f1, 2.0);
        if (isZero(f1)) {
            delete &f1;
            return f2;
        }
        if (isZero(f2)) {
            delete &f2;
            return f1;
        }
        doublereal c = f1.isProportional(f2);
        if (c != 0) {
            if (c == -1.0) 
                return *(new Const1(0.0));
            else {
                return newTimesConstFunction(f1, c + 1.0);
            }
        }
        return *(new Sum1(f1, f2));
    }

    Func1& newDiffFunction(Func1& f1, Func1& f2) {
        if (isZero(f2)) {
            delete &f2; return f1;
        }
        if (f1.isIdentical(f2)) {
            delete &f1; delete &f2; 
            return *(new Const1(0.0));
        }
        doublereal c = f1.isProportional(f2);
        if (c != 0.0) {
            if (c == 1.0) return *(new Const1(0.0));
            else return newTimesConstFunction(f1, 1.0 - c);
        }
        return *(new Diff1(f1, f2));
    }

    Func1& newProdFunction(Func1& f1, Func1& f2) {
        if (isOne(f1)) {
            delete &f1; return f2;
        }
        if (isOne(f2)) {
            delete &f2; return f1;
        }
        if (isZero(f1) || isZero(f2)) {
            delete &f1; delete &f2;
            return *(new Const1(0.0));
        }
        if (isConstant(f1) && isConstant(f2)) {
            doublereal c1c2 = f1.c() * f2.c();
            delete &f1; delete &f2;
            return *(new Const1(c1c2));
        }
        if (isConstant(f1)) {
            doublereal c = f1.c();
            delete &f1;
            return newTimesConstFunction(f2, c);
        }
        if (isConstant(f2)) {
            doublereal c = f2.c();
            delete &f2;
            return newTimesConstFunction(f1, c);
        }
            
        if (isPow(f1) && isPow(f2)) {
            Func1& p = *(new Pow1(f1.c() + f2.c()));
            delete &f1; delete &f2;
            return p;
        }

        if (isExp(f1) && isExp(f2)) {
            Func1& p = *(new Exp1(f1.c() + f2.c()));
            delete &f1; delete &f2;
            return p;
        }

        bool tc1 = isTimesConst(f1);
        bool tc2 = isTimesConst(f2);

        if (tc1 || tc2) {
            doublereal c1 = 1.0, c2 = 1.0;
            Func1 *ff1 = 0, *ff2 = 0;
            if (tc1) {
                c1 = f1.c();
                ff1 = &f1.func1_dup();
                delete &f1;
            }
            else ff1 = &f1;
            if (tc2) {
                c2 = f2.c();
                ff2 = &f2.func1_dup();
                delete &f2;
            }
            else ff2 = &f2;
            Func1& p = newProdFunction(*ff1, *ff2);

            if (c1*c2 != 1.0) {
                return newTimesConstFunction(p, c1*c2);
            }
            else
                return p;
        }
        else
            return *(new Product1(f1, f2));
    }

    Func1& newRatioFunction(Func1& f1, Func1& f2) {
        if (isOne(f2)) return f1;
        if (isZero(f1)) return *(new Const1(0.0));
        if (f1.isIdentical(f2)) {
            delete &f1; delete &f2;
            return *(new Const1(1.0));
        }
        if (f1.ID() == PowFuncType && f2.ID() == PowFuncType) {
            return *(new Pow1(f1.c() - f2.c()));
        }
        if (f1.ID() == ExpFuncType && f2.ID() == ExpFuncType) {
            return *(new Exp1(f1.c() - f2.c()));
        }
        return *(new Ratio1(f1, f2));
    }
    
    Func1& newCompositeFunction(Func1& f1, Func1& f2) {
        //#ifdef DEBUG_FUNC
        //cout << "creating new composite function." << endl;
        //cout << "f1 = " << f1.write("x") << "  " << f1.ID() << endl;
        //cout << "f2 = " << f2.write("x") << "  " << f2.ID() << endl;
        //#endif
        if (isZero(f1)) {
            delete &f1; delete &f2;
            return *(new Const1(0.0));
        }
        if (isConstant(f1)) {
            delete &f2;
            return f1;
        }
        if (isPow(f1) && f1.c() == 1.0) {
            delete &f1;
            return f2;
        }
        if (isPow(f1) && f1.c() == 0.0) {
            delete &f1;
            delete &f2;
            return *(new Const1(1.0));
        }
        if (isPow(f1) && isPow(f2)) {
            doublereal c1c2 = f1.c() * f2.c(); 
            delete &f1;
            delete &f2;
            return *(new Pow1(c1c2));
        }
        return *(new Composite1(f1, f2));
    }

    Func1& newTimesConstFunction(Func1& f, doublereal c) {
        if (c == 0.0) {
            delete &f; 
            return *(new Const1(0.0));
        }
        if (c == 1.0) {
            return f;
        }
        if (f.ID() == TimesConstantFuncType) {
            f.setC(f.c() * c);
            return f;
        }
        return *(new TimesConstant1(f, c));
    }

    Func1& newPlusConstFunction(Func1& f, doublereal c) {
        if (c == 0.0) {
            return f;
        }
        if (isConstant(f)) {
            doublereal cc = f.c() + c;
            delete &f;
            return *(new Const1(cc));
        }
        if (f.ID() == PlusConstantFuncType) {
            f.setC(f.c() + c);
            return f;
        }
        return *(new PlusConstant1(f, c));
    }

}
