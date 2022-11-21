/**
 * @file factory.hpp
 * @author M. Puetz
 * @brief Generic class for object factories.
 * @date 2022-09-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <map>
#include <memory>
#include <functional>


/**
 * @brief Generic class for object factories.
 * 
 * @tparam BaseType Base type of the objects to be created.
 * @tparam ArgTypes Argument types of the constructor called during creation of objects.
 */
template<class BaseType, class... ArgTypes>
class Factory {

private:

        /**
         * @brief Key that is used for an element that is 'None' corresponding to a null-pointer.
         * 
         */
        static const std::string& noneKey_()
        {
            static std::string noneKey = "None";
            return noneKey;
        }

        /**
         * @brief Vector of keys (typically the registered type names).
         * 
         */
        static std::vector<std::string>& keys_()
        {
            static std::vector<std::string> keys = {};
            return keys;
        }

        /**
         * @brief Map of type names and the respective functions for the creation of objects.
         * 
         */
        static auto& map_()
        {
            static std::map<std::string, std::function<BaseType* (ArgTypes...)> > map = {};
            return map;
        }


public:

        /**
         * @brief Create object of the type corresponding to the given string and return unique pointer.
         * 
         * @param key Key of the type to be instantiated corresponding to `map_`.
         * @param args Constructor arguments.
         * @return std::unique_ptr<BaseType> Unique pointer to created object.
         */
        static std::unique_ptr<BaseType> makeUnique(const std::string& key, ArgTypes... args) {

            if (key == noneKey_()) {
                return std::unique_ptr<BaseType>(nullptr);
            }

            try {
                return std::unique_ptr<BaseType>(map_()[key](args...));
            }
            catch (std::bad_function_call &e) {
                std::string errMsg = "The type '" + key 
                    + "' does probably not exist. Original exception: " + e.what();
                throw std::runtime_error(errMsg);
            }
        }

        /**
         * @brief Create object of the type corresponding to the given string and return shared pointer.
         * 
         * @param key Key of the type to be instantiated corresponding to `map_`.
         * @param args Constructor arguments.
         * @return std::shared_ptr<BaseType> Shared pointer to created object.
         */
        static std::shared_ptr<BaseType> makeShared(const std::string& key, ArgTypes... args) {

            if (key == noneKey_()) {
                return std::shared_ptr<BaseType>(nullptr);
            }

            try {
                return std::shared_ptr<BaseType>(map_()[key](args...));
            }
            catch (std::bad_function_call &e) {
                std::string errMsg = "The type '" + key 
                    + "' does probably not exist. Original exception: " + e.what();
                throw std::runtime_error(errMsg);
            }
        }

        /**
         * @brief Register type derived from `BaseType`.
         * 
         * @param key Key corresponding to the registered type, which is recommended to be the type name.
         * @param creatorFunc Function returning a poiner to new object of registered type.
         * @return true If type has been registered successfully.
         * @return false If registration of the type failed because the key already exists.
         */
        static bool registerType(const std::string &key, 
                std::function<BaseType* (ArgTypes...)> creatorFunc) {
            
            bool keyExists = static_cast<bool>(map_().count(key));
            if (!keyExists) {
                map_()[key] = creatorFunc;
                keys_().push_back(key);
            }

            return !keyExists;
        }

        /**
         * @brief Register type derived from `BaseType`.
         * 
         * @tparam SubType The type derived from `BaseType` that is supposed to be registered
         * @param key Key corresponding to the registered type, which is recommended to be the type name.
         * @return true If type has been registered successfully.
         * @return false If registration of the type failed because the key already exists.
         */
        template<class SubType>
        static bool registerType(const std::string &key){

            std::function<BaseType* (ArgTypes...)> 
                creatorFunc = [](ArgTypes... args){return new SubType(args...);};

            return registerType(key, creatorFunc);
        }

        /**
         * @brief Unregister type.
         * 
         * @param key Key of the type to be unregistered.
         * @return true If key was found.
         * @return false If key was not found.
         */
        static bool unregisterType(const std::string &key) {

            bool foundKey = static_cast<bool>(!map_().erase(key));

            if (foundKey) {
                keys_().clear();
                for (auto& item : map_) {
                    keys_().push_back(item.second);
                }
            }

            return foundKey;
        }

        /**
         * @brief Get constant reference to the vector of keys.
         * 
         * @return const std::vector<std::string>& Vector of keys.
         */
        static const std::vector<std::string>& keys() {
            return keys_();
        }

        static const std::string& noneKey() {
            return noneKey_();
        }

};

#define REGISTER_TYPE(FACTORY_NAME, TYPE_NAME) \
        FACTORY_NAME::registerType<TYPE_NAME>(#TYPE_NAME)

#endif // FACTORY_HPP