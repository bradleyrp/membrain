res_collection = pickle.load(open('/home/rpb/worker/restemp','r'))
plt.imshow(mean([i[0] for i in res_collection],axis=0));plt.show()
plt.hist([i for j in mean([i[0] for i in res_collection],axis=0) for i in j]);plt.show()
pickle.dump(res_collection,open('/home/rpb/v550shortrescol','w'))
