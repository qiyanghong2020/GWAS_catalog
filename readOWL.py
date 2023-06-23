from owlready2 import *

# 加载OWL文件
onto = get_ontology("efo.owl").load()

# 获取指定类的子类
parent_class = onto.get_iri_for_abbrev("ParentClass")
subclasses = list(parent_class.subclasses())

# 打印子类
for subclass in subclasses:
    print("Subclass:", subclass)
    

# 访问本体中的类
for cls in onto.classes():
    print("Class:", cls)

# 访问本体中的属性
for prop in onto.object_properties():
    print("Object Property:", prop)

for prop in onto.data_properties():
    print("Data Property:", prop)

# 访问本体中的实例
for ind in onto.individuals():
    print("Individual:", ind)
