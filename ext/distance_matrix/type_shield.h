/* 
 * File:   type_shield.h
 * Author: jwoods
 *
 * This is a template for a special basic type which can "shield" another type.
 * It behaves as its first template argument, but that shields another hidden
 * value which can be accessed by reference using get().
 *
 * Be aware that subtraction and addition operations are not implemented, as
 * the shield's behavior would be undefined. Nevertheless, comparison operators
 * are available; these refer only to the val, not the shield.
 *
 * Created on July 23, 2009, 11:46 AM
 */

#ifndef _TYPE_SHIELD_H
#define	_TYPE_SHIELD_H

template <typename T, typename S>
class type_shield {
public:
    type_shield(const T& t, const S& s) : val(t), shield(s) { }
    type_shield(const type_shield<T,S>& orig) : val(orig.val), shield(orig.shield) { }
    virtual ~type_shield() { }
    type_shield<T,S>& operator=(const type_shield<T,S>& orig) { val = orig.val; shield = orig.shield; }
    

    // Comparison operators -- only compare val, not shield!
    bool operator==(const type_shield<T,S>& rhs) const {
        return (*this == rhs.val); // defer to the == T& rhs operator.
    }
    bool operator==(const T& rhs) const {
        return (val == rhs);
    }
    bool operator!=(const type_shield<T,S>& rhs) const {
        return (*this != rhs.val); // defer
    }
    bool operator!=(const T& rhs) const {
        return (val != rhs);
    }
    bool operator<=(const type_shield<T,S>& rhs) const {
        return (*this <= rhs.val);
    }
    bool operator<=(const T& rhs) const {
        return (val <= rhs);
    }
    bool operator>=(const type_shield<T,S>& rhs) const {
        return (*this >= rhs.val);
    }
    bool operator>=(const T& rhs) const {
        return (val >= rhs);
    }
    bool operator<(const type_shield<T,S>& rhs) const {
        return (*this < rhs.val);
    }
    bool operator<(const T& rhs) const {
        return (val < rhs);
    }
    bool operator>(const type_shield<T,S>& rhs) const {
        return (*this > rhs.val);
    }
    bool operator>(const T& rhs) const {
        return (val > rhs);
    }

    // Operators that return a T instead of a type_shield
    T operator*(const T& rhs) const {
        return val * rhs;
    }
    T operator/(const T& rhs) const {
        return val / rhs;
    }


    // Accessors for shielded value
    S& get() { return shield; }
    const S& get() const { return shield; }

    // Also have it behave somewhat like a pair.
    T& first() { return val; }
    const T& first() const { return val; }
    S& second() { return shield; }
    const S& second() const { return shield; }
    
protected:
    T val;
    S shield;
};

//#include "type_shield.cpp"

#endif	/* _TYPE_SHIELD_H */

